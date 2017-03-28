#include "file.reading.hh"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "other/DIE.hh"
#include "other/PP.hh"
#include "other/utils.hh"

using std:: ifstream;
using std:: string;
using utils:: ssize;

#include "fwd/src/file.reading.hh"

namespace file_reading {

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
    int pos = 0;
    while(1) {
        auto next_delim = header_line.find(delimiter, pos);
        PP(header_line.substr(pos, next_delim-pos));
        if(next_delim== string:: npos)
            break;
        pos=next_delim+1;
    }
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
