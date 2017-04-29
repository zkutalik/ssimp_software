#include "compression.hh"

#include <iostream>
#include <fstream>

#include <gzstream.h>

#include "other/DIE.hh"

using std:: ofstream;
using std:: ios_base;
using std:: string;

namespace compression {
    void make_compressed_vcf_file   (std:: string compressed_out_file_name, std:: string vcf_filename) {
        ofstream binary_output  (compressed_out_file_name
                                ,   ios_base::out
                                 |  ios_base::binary
                                 |  ios_base::trunc
                                );
        (!!binary_output) || DIE("Couldn't open [" << compressed_out_file_name.c_str()
                << "] for writing.");

        gz:: igzstream f(vcf_filename.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << vcf_filename << ']');
        string current_line;

        //int simple_line_number = 0;
    }
}
