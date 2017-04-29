#include "compression.hh"

#include <iostream>
#include <fstream>

#include <gzstream.h>

#include "other/DIE.hh"
#include "other/range_action.hh"
#include "other/range_view.hh"
#include "other/PP.hh"

namespace action = range:: action;
namespace view   = range:: view  ;

using std:: ofstream;
using std:: ios_base;
using std:: string;

namespace range {
    // put this into the 'range' namespace, to help begin() find it via ADL
    struct read_file_as_a_range_t {
        gz:: igzstream      &   m_f;
        string                  current_line;

        read_file_as_a_range_t(gz:: igzstream & f) : m_f(f) {}

        bool            empty()     const   { return !m_f; }
        void            advance()           { getline(m_f, current_line); }
        const string&   front_ref() const   { return current_line; }
        string          front_val() const   { return current_line; }
    };
}
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

        range:: read_file_as_a_range_t vcf_input_range{f};
        auto z = zip( range::ints(10000000), vcf_input_range | view:: ref_wraps );

        //for(string const & s : vcf_input_range) { PP(s.substr(0,20)); }
        //z |action:: foreach| [&](auto) { };
        for(auto const & x : z) {
            PP( std:: get<0>(x)
              , std:: get<1>(x).get().substr(0,20) );
        }

    }
}
