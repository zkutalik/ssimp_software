#include "file.reading.hh"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "other/DIE.hh"
#include "other/PP.hh"

using std:: ifstream;
using std:: string;

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
    decide_delimiter(header_line);
}

FWD(file_reading)
static
char   decide_delimiter( string      const & header_line ) {
    // Which is comma, tab, or space, are most common here?

    size_t commas = std::count(header_line.begin(), header_line.end(), ',');
    size_t tabs   = std::count(header_line.begin(), header_line.end(), '\t');
    size_t spaces = std::count(header_line.begin(), header_line.end(), ' ');
    PP(commas, tabs, spaces);
}

} // namespace file_reading
