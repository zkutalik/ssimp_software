#include "file.reading.hh"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cassert>

#include "other/DIE.hh"
#include "other/PP.hh"
#include "other/utils.hh"

using std:: ifstream;
using std:: string;
using std:: vector;
using utils:: ssize;
using utils:: operator<<; // to print vectors

#include "fwd/src/file.reading.hh"

namespace file_reading {

struct header_details {
    struct offset_and_name {
        int                  m_offset;
        string               m_name;

        offset_and_name() : m_offset(-1)
        {}
        offset_and_name(int offset, string name) : m_offset(offset), m_name(name)
        {}
        offset_and_name(offset_and_name const &) = default;
        offset_and_name(offset_and_name      &&) = default;

        void operator=(offset_and_name const & from) {
            if(m_offset != -1) {
                WARNING("Two fields with a similar name, using [" << m_name << "] and ignoring [" << from.m_name << "].");
                return;
            }
            assert(m_offset == -1);
            m_offset = from.m_offset;
            m_name   = from.m_name;
        }
        void operator=(offset_and_name && from) {
            return *this = from;
        }
    };
    offset_and_name    SNPname;
    offset_and_name    chromosome;
    offset_and_name    position;
    offset_and_name    allele_ref;
    offset_and_name    allele_alt;
    offset_and_name    qual;
    offset_and_name    filter;
    offset_and_name    info;
    offset_and_name    format;
    vector<offset_and_name> unaccounted;
};

struct PlainVCFfile : public Genotypes_I {
};

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name) {
    // I really should detect the file-type.
    // But for now, we'll just assume a plain (non-gzipped) vcf file.
    return read_in_a_raw_ref_file_as_VCF(file_name);
}

FWD(file_reading)
static
GenotypeFileHandle      read_in_a_raw_ref_file_as_VCF(std:: string file_name) {
    PP(file_name);
    ifstream f(file_name);
    string line;
    while(getline(f, line)) {
        if(line.at(0) == '#' && line.at(1) == '#')
            continue; // skip these ## lines
        parse_header(line);
        break;
    }
    return {};
}

FWD(file_reading)
static
void   parse_header( string      const & header_line ) {
    PP(header_line);
    char delimiter = decide_delimiter(header_line);
    auto field_names = tokenize(header_line, delimiter);
    PP(field_names.size());
    PP(field_names);

    header_details hd;

    int field_counter = -1;
    for(auto & one_field_name : field_names) {
        ++field_counter;
        // go through each field name in turn and try to account for it
        if(false) {}
        else if(is_in_this_list(one_field_name, {"ID"})) {
            hd.SNPname = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"#CHROM","chr"})) {
            hd.chromosome = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"POS"})) {
            hd.position = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"REF"})) {
            hd.allele_ref = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"ALT"})) {
            hd.allele_alt = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"QUAL"})) {
            hd.qual = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"FILTER"})) {
            hd.filter = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"INFO"})) {
            hd.info = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"FORMAT"})) {
            hd.format = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else {
            hd.unaccounted.push_back( header_details:: offset_and_name(field_counter, one_field_name) );
        }
    }
    auto checker = [&](auto & offset_and_name) {
        if(offset_and_name.m_offset == -1) {
            assert(offset_and_name.m_name == "");
            return;
        }
        PP(offset_and_name.m_name);
        PP(field_names.at(offset_and_name.m_offset));
        assert(offset_and_name.m_name == field_names.at(offset_and_name.m_offset));
    };
    checker(hd.SNPname);
    checker(hd.chromosome);
    checker(hd.position);
    checker(hd.allele_ref);
    checker(hd.allele_alt);
    checker(hd.qual);
    checker(hd.filter);
    checker(hd.info);
    checker(hd.format);
    PP(hd.unaccounted.size());
}

FWD(file_reading)
static bool is_in_this_list(string const & s, std:: initializer_list<char const *> candidates) {
    for(auto cand : candidates) {
        if(s==cand)
            return true;
    }
    return false;
}

FWD(file_reading)
static
vector<string>   tokenize( string      const & line
               , char                delimiter
        ) {
    int pos = 0;
    vector<string> fields;
    while(1) {
        auto next_delim = line.find(delimiter, pos);
        fields.push_back( line.substr(pos, next_delim-pos) );
        if(next_delim== string:: npos)
            break;
        pos=next_delim+1;
    }
    return fields;
}

FWD(file_reading)
static
char   decide_delimiter( string      const & header_line ) {
    // Which is comma, tab, or space, are most common here?

    size_t commas = std::count(header_line.begin(), header_line.end(), ',');
    size_t tabs   = std::count(header_line.begin(), header_line.end(), '\t');
    size_t spaces = std::count(header_line.begin(), header_line.end(), ' ');
    char delimiter = '\0';
    if(commas >  0 && tabs == 0 && spaces == 0) delimiter = ',';
    if(commas == 0 && tabs >  0 && spaces == 0) delimiter = '\t';
    if(commas == 0 && tabs == 0 && spaces >  0) delimiter = ' ';
    if(delimiter == '\0') {
        DIE("Couldn't decide what delimiter to use."
                <<  " #commas:" << commas
                << ", #tabs:"     << tabs
                << ", #spaces:"     << spaces
           );
    }
    return delimiter;
}

} // namespace file_reading
