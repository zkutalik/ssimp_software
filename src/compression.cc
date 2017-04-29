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

using utils:: operator<<;

namespace range {
    // put this into the 'range' namespace, to help begin() find it via ADL
    struct read_file_as_a_range_t {
        gz:: igzstream      &   m_f;
        string                  current_line;

        read_file_as_a_range_t(gz:: igzstream & f) : m_f(f) {
            advance();
        }

        bool            empty()     const   { return !m_f; }
        void            advance()           { getline(m_f, current_line); }
        const string&   front_ref() const   { return current_line; }
        string          front_val() const   { return current_line; }
    };
}
namespace compression {

    struct GTcompressed_output_t {
        ofstream m_f;

        template<typename ...Ts>
        GTcompressed_output_t(Ts&& ...ts) : m_f( std::forward<Ts>(ts) ... ) {}

        operator bool() const {
            return !!m_f;
        }
    };

    void make_compressed_vcf_file   (std:: string compressed_out_file_name, std:: string vcf_filename) {
        GTcompressed_output_t binary_output (compressed_out_file_name
                                            ,   ios_base::out
                                               |ios_base::binary
                                               |ios_base::trunc
                                            );
        (  binary_output) || DIE("Couldn't open [" << compressed_out_file_name.c_str()
                << "] for writing.");

        gz:: igzstream f(vcf_filename.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << vcf_filename << ']');

        range:: read_file_as_a_range_t vcf_input_range{f};
        auto z = zip( range::ints<int64_t>(), vcf_input_range );

        constexpr int INFO_field_to_skip = 7;

        std::move(z) |action:: unzip_foreach|
        [&](auto i, std:: string const & s) {
            if(s.substr(0,2) == "##")
                return;
            auto fields = utils:: tokenize(s, '\t');
            //PP(fields);
            if(fields.at(0) == "#CHROM") {
                // #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
                assert( fields.at(INFO_field_to_skip) == "INFO" );
                assert((
                std::vector<string> {   fields.at(0)
                                    ,   fields.at(1)
                                    ,   fields.at(2)
                                    ,   fields.at(3)
                                    ,   fields.at(4)
                                    ,   fields.at(5)
                                    ,   fields.at(6)
                                    ,   fields.at(7)
                                    ,   fields.at(8)
                                    } ==
                std::vector<string> { "#CHROM" , "POS" , "ID" , "REF", "ALT", "QUAL", "FILTER" , "INFO", "FORMAT"
                                    }));
                return;
            }
            PP( i, s.substr(0,90) );
        };

    }
}
