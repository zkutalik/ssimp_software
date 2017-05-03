#include<iostream>
#include<sstream>

#include "gzstream.h"

#include "range/range_from.hh"
#include "range/range_view.hh"
#include "range/range_action.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "zlib-vector-of-char/zlib-vector.hh"

namespace from = range:: from;
namespace action = range:: action;
namespace view   = range:: view  ;

using std:: array;
using std:: vector;
using std:: string;
using std:: ofstream;
using std:: cout;
using std:: ios_base;
using std:: ostringstream;
using utils:: tokenize;
using utils:: stdget0;
using utils:: stdget1;
using utils:: print_type;
using utils:: nice_operator_shift_left;
using utils:: save_ostream_briefly;
using range:: range_from_begin_end;
using range:: view:: take;

enum class vcfGTz_codes : uint8_t {
            code_end_of_file        = 0
        ,   code_null_term_string   = 1 // null-terminated string - when everything else fails
};

struct vcfGTz_writer {
    ofstream m_f;
    decltype(m_f.tellp()) my_beginning_pos;

    vcfGTz_writer(string output_filename) : m_f {   output_filename
                                                ,   ios_base::out
                                                  | ios_base::binary
                                                  | ios_base::trunc
    }{
        my_beginning_pos = m_f.tellp();

        constexpr char magic_file_header[] = "vcfGTz.0.0.1";
        save_string0( magic_file_header );

        assert( sizeof(magic_file_header) == get_current_file_offset() );
    }

    void            save_this_line(vector<string>   const & fields) {
        ensure_there_are_no_nulls(fields);

        constexpr int HOW_MANY_FIELDS_TO_SAVE_DIRECTLY = 7;
        // store the first seven fields directly, so these
        // can be used for lookups later with having to go via
        // zlib
        // Then compress the remainder of the fields using zlib

        for(int f = 0; f < HOW_MANY_FIELDS_TO_SAVE_DIRECTLY; ++f) {
            save_smart_string0(fields.at(f));
        }

        auto just_last_fields = range:: range_from_begin_end (
                fields.begin() + HOW_MANY_FIELDS_TO_SAVE_DIRECTLY
                ,fields.end()
                );

        ostringstream detokenized_stream;
        for(auto & any_field : just_last_fields) {
            if(&any_field != &fields.front()) {
                detokenized_stream << '\t';
            }
            detokenized_stream << any_field;
        }
        string detokenized = detokenized_stream.str();

        zlib_vector:: vec_t as_a_vector( detokenized.begin(), detokenized.end() );
        zlib_vector:: vec_t compressed = zlib_vector:: deflate( as_a_vector );

        save_vector_of_char_with_leading_size( compressed );

        {
            //assert(fields == tokenize(detokenized,'\t'));
            assert(as_a_vector == zlib_vector:: inflate( compressed ));
        }
    }

    void save_smart_string0(    string const & s) {
        save_code(vcfGTz_codes:: code_null_term_string);
        save_string0(s.c_str());
    }

    void            save_string0(char const *s) {
        m_f << s;
        m_f << '\0';
    }
    void            save_vector_of_char_with_leading_size( zlib_vector:: vec_t const & compressed ) {
        uint32_t sz = compressed.size();
        save_uint32_t(sz);

        static_assert(std::is_same< unsigned char const * , decltype(compressed.data()) >{} ,"");
        m_f.write   (   reinterpret_cast<char const *>(compressed.data())
                    ,   compressed.size()
                    );
    }

    struct bit_conversions {
        template<typename T>
    static
    array<uint8_t,4> convert_to_four_bytes(T u) {
        static_assert(std:: is_same<T, uint32_t>{}, "");

        //cout << '\n'; PP(u);
        {
            save_ostream_briefly sob{cout};
            cout << std::hex << u << '\n';
        }
        array<uint8_t,4> arr;
        static_assert(4==utils::ssize(arr) ,"");
        for(int byte = 0; byte < utils::ssize(arr); ++byte) {
            auto by = (u >> (byte*8)) & 255;
            arr.at(byte) = by;
        }
        return arr;
    }
    static
    uint32_t convert_from_four_bytes(array<uint8_t,4> arr) {
        static_assert(arr.size() == 4 ,""); // to return a uint32_t
        uint32_t u = 0;
        for(int byte = arr.size(); byte > 0; --byte) {
            u <<= 8;
            u += arr.at(byte-1);
        }
        return u;
    }
    };

    void            save_uint32_t(uint32_t sz) {
        auto arr4bytes = bit_conversions:: convert_to_four_bytes(sz);
        assert(sz == bit_conversions:: convert_from_four_bytes(arr4bytes));
        for(uint8_t u8: arr4bytes) {
            m_f << (char) u8;
        }
    }

    void            save_code(vcfGTz_codes code) {
        m_f << (char)(uint8_t)code;
    }

    int64_t         get_current_file_offset() {
        // not every platform uses anything as big as int64_t,
        // but I'm going to force it, as I want something
        // standard that I can store in the file
        assert(m_f);
        auto off = m_f.tellp() - my_beginning_pos;
        assert(off > 0);
        return off;
    }
    void            ensure_there_are_no_nulls(vector<string>   const & fields) {
        for(auto & f : fields) {
            for(char const & c : f) {
                assert( c!='\0' );
            }
        }
    }
};

int main(int argc, char **argv) {
    // Compress the input vcf file. The following info is discarded from the vcf file:
    //  - The 'preamble'. Everything before before the header ("#CHROM\tPOS"...) is ignored
    //      * (We do keep the header though)
    //  - The INFO field is ignored - too long and not very interesting
    //  - only the 'GT' data is kept, i.e. the other four in 'GT:PL:DP:GP:PS' are ignored is present
    // Everything else is faithfully compressed (via 'zlib') line-by-line and stored.
    // An index is appended to allow very fast access to the input in three orders:
    //  - original file order
    //  - CHROM:POS order
    //  - ID order
    argc == 3 || DIE(argv[0] << ": requires 2 args, input (.vcf/.vcfgz) and output");

    string arg_input_filename = argv[1];
    string arg_output_filename = argv[2];

    gz:: igzstream f(arg_input_filename.c_str());
    f.fail() && DIE("Can't find file [" << arg_input_filename << ']');

    auto r = zip(from::ifstream( f ), range:: ints()) ;

    while(!r.empty() && std::get<0>(range:: front_val(r)).substr(0,2) == "##")
        advance(r);

    int64_t SNP_counter = -2; // because -1 is the 'special' SNP for the header

    vcfGTz_writer writer{arg_output_filename};

    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        auto && s = x |stdget0;

        ++SNP_counter;
        auto counter_at_this_SNP = SNP_counter;

        auto fields = tokenize(s, '\t');
        PP(counter_at_this_SNP, fields.size());

        if(SNP_counter == -1) {
            auto first_9_fields =
            range_from_begin_end(fields) |take|9 |action::collect;
            first_9_fields == tokenize("#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT",',') || DIE("Malformed header line?" << fields);
        }

        // Copy the two interesting fields before deleting them
        // ** all lines, not just the header **

        constexpr size_t index_of_field_INFO = 7; // MUST be three consecutive values
        constexpr size_t index_of_field_FORMAT = 8; // MUST be three consecutive values
        constexpr size_t index_of_field_therest = 9; // MUST be three consecutive values
        constexpr size_t index_of_field_therest_after_deletion = 7;

        auto INFO = fields.at(index_of_field_INFO);
        auto FORMAT = fields.at(index_of_field_FORMAT);
        fields.erase( fields.begin() + index_of_field_INFO, fields.begin() + index_of_field_therest );

        if(SNP_counter >= 0) {
            // Next, delete the non-'GT' fields from the individual data (and from the header, for completeness)
            auto subfields_that_are_GT =
                from :: vector(tokenize(FORMAT,':'))
                |view::map| [](auto && s){ return s == "GT"; }
                |view:: which
                |action:: collect
            ;
            subfields_that_are_GT.size() == 1 || DIE(
                    "Was expecting 'GT' exactly once in the format field. ["
                    << arg_input_filename << ':' << line_no
                    << "] ["
                    << FORMAT
                    << "]"
                    );

            auto subfield_index_for_GT = subfields_that_are_GT.at(0);
            (void)subfield_index_for_GT;
            auto r = range_from_begin_end(  fields.begin() + index_of_field_therest_after_deletion
                                         ,  fields.end()
                                        );
            while(!r.empty()) {
                auto & x = front_ref(r);
                auto t = tokenize(x,':');
                (subfield_index_for_GT < t.size()) || DIE("FORMAT problem at "
                    << '[' << arg_input_filename << ':' << line_no << ']'
                    << " [" << FORMAT << "] [" << x << "]");
                x=
                t.at(subfield_index_for_GT);

                r.advance();
            }

        }

        /* Now we've finally got a version of the line worth saving
         * Remember that this includes the header line
         */
        //PP(fields.size(), fields);
        writer.save_this_line(fields);
    }
}
