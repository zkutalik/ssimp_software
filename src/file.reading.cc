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

#define LOOKUP(hd, fieldname, vec) lookup(hd . fieldname, vec, #fieldname)

namespace file_reading {
    SNPiterator &       SNPiterator:: operator++()        {
        assert(m_line_number < m_gfh->number_of_snps());
        ++m_line_number;
        return *this;
    }
    bool                SNPiterator:: operator==(SNPiterator const & other) {
        assert(m_gfh == other.m_gfh);
        return m_line_number == other.m_line_number;
    }
    bool                SNPiterator:: operator!=(SNPiterator const & other) {
        return !(*this == other);
    }
    bool                SNPiterator:: operator< (SNPiterator const & other) {
        assert(m_gfh == other.m_gfh);
        return m_line_number < other.m_line_number;
    }
    bool                SNPiterator:: operator>=(SNPiterator const & other) { return !(*this < other); }
    int                 SNPiterator:: operator- (SNPiterator const & other) {
        assert(m_gfh == other.m_gfh);
        return m_line_number - other.m_line_number;
    }
    chrpos              SNPiterator:: operator* () const {
        return get_chrpos();
    }
    SNPiterator &       SNPiterator:: operator+=(long int            ran)   {
        m_line_number += ran;
        return *this;
    }

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

    char                     m_delimiter;
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

} // namespace file_reading

#include "fwd/src/file.reading.hh"

namespace file_reading {

struct OneLineSummary {
    ifstream:: pos_type      m_tellg_of_line_start;
    string                   m_SNPname;
    int                      m_chromosome;
    int                      m_position;
    string                   m_allele_alt;
    string                   m_allele_ref;
};
struct GwasLineSummary {
    string                   m_SNPname;
    //chrpos                   m_chrpos;
    string                   m_allele_alt;
    string                   m_allele_ref;
    double                   m_z;
};

struct PlainVCFfile : public file_reading:: Genotypes_I
{
    vector<OneLineSummary> m_each_SNP_and_its_offset;
    string                 m_underlying_file_name;

    virtual int         number_of_snps() const {
        return m_each_SNP_and_its_offset.size();
    }
    virtual chrpos      get_chrpos         (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return chrpos{ols.m_chromosome, ols.m_position};
    }
    virtual std::string get_SNPname        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_SNPname;
    }
    virtual std::string get_allele_ref        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_allele_ref;
    }
    virtual std::string get_allele_alt        (int i)     const {
        OneLineSummary const & ols = get_ols(i);
        return ols.m_allele_alt;
    }

    OneLineSummary  get_ols         (int i)     const {
        assert(i>=0);
        assert(i<number_of_snps());
        return m_each_SNP_and_its_offset.at(i);
    }

};

GenotypeFileHandle      read_in_a_raw_ref_file(std:: string file_name) {
    // I really should detect the file-type.
    // But for now, we'll just assume a plain (non-gzipped) vcf file.
    return read_in_a_raw_ref_file_as_VCF(file_name);
}

GwasFileHandle      read_in_a_gwas_file(std:: string file_name) {
    // I really should detect the file-type.
    // But for now, we'll just assume a plain (non-gzipped) vcf file.
    return read_in_a_gwas_file_simple(file_name);
}


FWD(file_reading)
static
GenotypeFileHandle      read_in_a_raw_ref_file_as_VCF(std:: string file_name) {
    PP(file_name);
    ifstream f(file_name);
    string current_line;

    header_details hd;

    // Skip past the '##' lines
    while(getline(f, current_line)) {
        if(current_line.at(0) == '#' && current_line.at(1) == '#')
            continue; // skip these ## lines
        break;
    }

    // current_line is now the first line, i.e. the header
    hd = parse_header(current_line);


    auto p = std:: make_shared<PlainVCFfile>();
    p->m_underlying_file_name = file_name;

    while(1) {
        OneLineSummary ols;
        ols.m_tellg_of_line_start = f.tellg();
        getline(f, current_line);
        if(!f) {
            f.eof() || DIE("Error before reaching eof() in this file");
            break;
        }
        auto all_split_up = tokenize(current_line, hd.m_delimiter);
        try {
            ols.m_SNPname    =                           LOOKUP(hd, SNPname, all_split_up);
            ols.m_chromosome = utils:: lexical_cast<int>( LOOKUP(hd, chromosome, all_split_up) );
            ols.m_position   = utils:: lexical_cast<int>( LOOKUP(hd, position, all_split_up) );
            ols.m_allele_alt =                           LOOKUP(hd, allele_alt, all_split_up);
            ols.m_allele_ref =                           LOOKUP(hd, allele_ref, all_split_up);

            p->m_each_SNP_and_its_offset.push_back(ols);
        } catch (std:: invalid_argument &e) {
            WARNING( "Ignoring this line, problem with the chromosome and/or position. ["
                    << "SNPname:" << LOOKUP(hd, SNPname   , all_split_up)
                    << " chr:" << LOOKUP(hd, chromosome, all_split_up)
                    << " pos:" << LOOKUP(hd, position  , all_split_up)
                    << "]" );
        }
    };

    sort( p->m_each_SNP_and_its_offset.begin()
        , p->m_each_SNP_and_its_offset.end()
        , [](auto &l, auto &r) -> bool {
            if(l.m_chromosome < r.m_chromosome) return true;
            if(l.m_chromosome > r.m_chromosome) return false;
            if(l.m_position   < r.m_position  ) return true;
            if(l.m_position   > r.m_position  ) return false;
            return false;
        }
        );

    return p;
}

FWD(file_reading)
static
string          lookup( header_details:: offset_and_name const &on
                      , vector<string> const & all_split_up
                      , const char *fieldname
                      ) {
    int offset = on.m_offset;
    (offset >= 0 && offset < ssize(all_split_up)) || DIE("Cannot find ["<<fieldname<<"] field");

    return all_split_up.at(offset);
}

FWD(file_reading)
static
file_reading::
header_details   parse_header( string      const & header_line ) {
    header_details hd;

    hd.m_delimiter = decide_delimiter(header_line);
    auto header_names = tokenize(header_line, hd.m_delimiter);

    int field_counter = -1;
    for(auto & one_field_name : header_names) {
        ++field_counter;
        // go through each field name in turn and try to account for it
        if(false) {}
        else if(is_in_this_list(one_field_name, {"ID","rnpid"})) {
            hd.SNPname = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"#CHROM","chr"})) {
            hd.chromosome = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"POS"})) {
            hd.position = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"REF","a1"})) {
            hd.allele_ref = header_details:: offset_and_name(field_counter, one_field_name);
        }
        else if(is_in_this_list(one_field_name, {"ALT","a2"})) {
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
    return hd;
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

struct SimpleGwasFile : public file_reading:: Effects_I
{
    vector<GwasLineSummary> m_each_SNP_and_its_z;
    string                 m_underlying_file_name;

    virtual int         number_of_snps() const {
        return m_each_SNP_and_its_z.size();
    }
    virtual std::string get_SNPname        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_SNPname;
    }
    virtual std::string get_allele_ref        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_allele_ref;
    }
    virtual std::string get_allele_alt        (int i)     const {
        auto const & ols = get_gls(i);
        return ols.m_allele_alt;
    }

    GwasLineSummary  get_gls         (int i)     const {
        assert(i>=0);
        assert(i<number_of_snps());
        return m_each_SNP_and_its_z.at(i);
    }

};

FWD(file_reading)
static
GwasFileHandle      read_in_a_gwas_file_simple(std:: string file_name) {
    PP(file_name);
    ifstream f(file_name);
    string current_line;

    header_details hd;

    getline(f, current_line);

    // current_line is now the first line, i.e. the header
    hd = parse_header(current_line);

    auto p = std:: make_shared<SimpleGwasFile>();
    p->m_underlying_file_name = file_name;

    while(1) {
        GwasLineSummary gls;
        getline(f, current_line);
        if(!f) {
            f.eof() || DIE("Error before reaching eof() in this file");
            break;
        }
        auto all_split_up = tokenize(current_line, hd.m_delimiter);
        try {
            gls.m_SNPname    =                           LOOKUP(hd, SNPname, all_split_up);
            //gls.m_chromosome = utils:: lexical_cast<int>( LOOKUP(hd, chromosome, all_split_up) );
            //gls.m_position   = utils:: lexical_cast<int>( LOOKUP(hd, position, all_split_up) );
            gls.m_allele_alt =                           LOOKUP(hd, allele_alt, all_split_up);
            gls.m_allele_ref =                           LOOKUP(hd, allele_ref, all_split_up);

            p->m_each_SNP_and_its_z.push_back(gls);
        } catch (std:: invalid_argument &e) {
            WARNING( "Ignoring this line, problem with the chromosome and/or position. ["
                    << "SNPname:" << LOOKUP(hd, SNPname   , all_split_up)
                    //<< " chr:" << LOOKUP(hd, chromosome, all_split_up)
                    //<< " pos:" << LOOKUP(hd, position  , all_split_up)
                    << "]" );
        }
    };

    /* Note: Sort them by position here */

    return p;
}

} // namespace file_reading
