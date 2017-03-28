#include "file.reading.hh"

#include <iostream>

#include "other/DIE.hh"
#include "other/PP.hh"

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
GenotypeFileHandle      read_in_a_raw_ref_file_as_VCF(std:: string file_name) {
    PP(file_name);
    return {};
}


} // namespace file_reading
