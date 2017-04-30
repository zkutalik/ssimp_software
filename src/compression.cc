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
using std:: vector;
using std:: cout;
using std:: endl;

using utils:: operator<<;
using utils:: stdget0;
using utils:: stdget1;

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

    static
    auto build_an_efficient_dictionary_from_a_vector_of_strings(vector<string> const & v) {
        (void) v;
    }

    enum class TypesOfCodedOutput : uint8_t {
            code_end_of_file        = 0
        ,   code_int                = 1
        ,   code_rs_int             = 2
        ,   code_period             = 3 // "." occurs a lot
        ,   code_PASS               = 4 // "PASS" occurs a lot
        ,   code_null_term_string   = 5 // null-terminated string - when everything else fails
        ,   code_A                  = 12
        ,   code_T                  = 13
        ,   code_G                  = 14
        ,   code_C                  = 15
    };
    struct GTcompressed_output_t {
        ofstream m_f;
        ofstream:: pos_type m_remember_the_begining_position;

        template<typename ...Ts>
        GTcompressed_output_t(Ts&& ...ts)   : m_f( std::forward<Ts>(ts) ... )
                                            , m_remember_the_begining_position( m_f.tellp() )
        {
            // Next two lines to ensure we really are at the start.
            m_f.seekp(0, ios_base::beg);
            assert( m_remember_the_begining_position == m_f.tellp() );
        }

        ~GTcompressed_output_t() { m_f.close(); }

        operator bool() const {
            return !!m_f;
        }

        void output_uint32(uint32_t x) {
            constexpr size_t w = sizeof(x);
            static_assert( w == 4 ,"");
            vector<uint8_t> bytes(w);
            for(int i=w-1; i>=0; --i) {
                bytes.at(i) = x % 256;
                x = x / 256;
            }
            assert(x==0);
            for(int i=0; i<(int)w; ++i) {
                m_f << bytes.at(i);
            }
        }
        void output_uint64(uint64_t x) {
            static_assert( 8 == sizeof(uint64_t) ,"");
            vector<uint8_t> bytes(8);
            for(int i=7; i>=0; --i) {
                bytes.at(i) = x % 256;
                x = x / 256;
            }
            assert(x==0);
            for(int i=0; i<8; ++i) {
                m_f << bytes.at(i);
            }
        }
        void output_string(string const & s) {
            for(auto c : s) {
                static_assert( std:: is_same<decltype(c), char>{} ,"");
                assert(c != '\0');
            }
            m_f << s;
            m_f << '\0';
        }
        void output_smart_string(string const & s) {
            for(auto c : s) {
                static_assert( std:: is_same<decltype(c), char>{} ,"");
                assert(c != '\0');
            }
            try {
                static_assert( sizeof(int) == sizeof(int32_t) ,"");
                int as_int = utils:: lexical_cast<int>(s);
                assert(as_int >= 0); // TODO: Could fail sometime
                output_code(TypesOfCodedOutput:: code_int);
                output_uint32(as_int);
                return;
            } catch (std:: invalid_argument & e) { }
            try {
                if(s.substr(0,2) == "rs") {
                    string s_skip_first_two = s.substr(2, string:: npos);
                    static_assert( sizeof(int) == sizeof(int32_t) ,"");
                    int as_int = utils:: lexical_cast<int>(s_skip_first_two);
                    assert(as_int >= 0); // TODO: Could fail sometime
                    output_code(TypesOfCodedOutput:: code_rs_int);
                    output_uint32(as_int);
                    return;
                }
            } catch (std:: invalid_argument & e) { }
            if(s == "."      ) { output_code(TypesOfCodedOutput:: code_period); return; }
            if(s == "PASS"   ) { output_code(TypesOfCodedOutput:: code_PASS  ); return; }
            if(s == "A"      ) { output_code(TypesOfCodedOutput:: code_A     ); return; }
            if(s == "T"      ) { output_code(TypesOfCodedOutput:: code_T     ); return; }
            if(s == "G"      ) { output_code(TypesOfCodedOutput:: code_G     ); return; }
            if(s == "C"      ) { output_code(TypesOfCodedOutput:: code_C     ); return; }
            // Otherwise, we just need to output a plain null-terminated string
            output_code(TypesOfCodedOutput:: code_null_term_string);
            output_string(s);
        }
        void output_code(TypesOfCodedOutput code) {
                m_f << (uint8_t)code;
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

        binary_output.output_string( "vcf.GT.compressed.0.0.1" ); // later, this will be the offset just past the end
        binary_output.output_uint64( 1234 ); // later, this will be the offset just past the end

        range:: read_file_as_a_range_t vcf_input_range{f};
        auto z = zip( range::ints<int64_t>(), vcf_input_range );

        constexpr int INFO_field_to_skip = 7;
        constexpr int FORMAT_field_to_skip = 8;
        constexpr int first_GT_field = 9;

        std::move(z) |action:: unzip_foreach|
        [&](auto i, std:: string const & s) {
            if(s.substr(0,2) == "##")
                return;
            auto fields = utils:: tokenize(s, '\t');
            //PP(fields);
            if(fields.at(0) == "#CHROM") {
                // #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
                assert( fields.at(INFO_field_to_skip) == "INFO" );
                assert( fields.at(FORMAT_field_to_skip) == "FORMAT" );
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
            (void)i;
            //PP( i, s.substr(0,90) );
            for(int j=0; j<9; ++j) {
                if(j==INFO_field_to_skip)
                    continue;
                if(j==FORMAT_field_to_skip)
                    continue;
                binary_output.output_smart_string( fields.at(j) );
            }

            /* Next, even though the FORMAT field isn't explicitly encoded, we
             * do need to find 'GT' in it
             */
            auto FORMAT = fields.at(FORMAT_field_to_skip);
            auto where_is_GT =
                view:: enumerate_vector(utils:: tokenize(FORMAT,':'))
                |view:: unzip_filter|
                [](auto /*index*/, auto fmt) {
                        return "GT" == fmt;
                }
                | action:: collect;
            where_is_GT.size() == 1 || DIE("The FORMAT field should have 'GT' exactly once [" << FORMAT << "]");

            int const offset_of_GT_within_FORMAT_field = where_is_GT.at(0) |stdget0;
            assert("GT" ==                              (where_is_GT.at(0) |stdget1));

            auto many_call_pairs_as_strings =
                range:: ints( first_GT_field, fields.size() )
                |view:: map|
                [&](auto index) {
                    return utils::tokenize(fields.at(index),':').at(offset_of_GT_within_FORMAT_field);
                }
                |action:: collect
                ;
            build_an_efficient_dictionary_from_a_vector_of_strings(many_call_pairs_as_strings);

        };

    }
}
