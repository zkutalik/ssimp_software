#include<iostream>
#include<sstream>

#include "gzstream.h"

#include "range/range_from.hh"
#include "range/range_view.hh"
#include "range/range_action.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "zlib-vector-of-char/zlib-vector.hh"

#include "vcfGTz_reader.hh"

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
using utils:: ELAPSED;
using range:: range_from_begin_end;
using range:: view:: take;
namespace bit_conversions = vcfGTz:: bit_conversions;

struct vcfGTz_writer {
    ofstream m_f;
    decltype(m_f.tellp()) my_beginning_pos;

    vcfGTz_writer(string output_filename) : m_f {   output_filename
                                                ,   ios_base::out
                                                  | ios_base::binary
                                                  | ios_base::trunc
    }{
        my_beginning_pos = m_f.tellp();

        constexpr char magic_file_header[] = "vcfGTz.0.0.2";
        save_string0( magic_file_header );

        assert( sizeof(magic_file_header) == get_current_file_offset() );
    }

    auto            start_a_new_block() {
        auto current_pos = m_f.tellp();
        uint64_t minus1 = ~ 0;
        save_uint64_t(minus1);
        return current_pos;
    }
    void            close_this_block( decltype(my_beginning_pos) pos_to_write_in_offset ) {
        auto current_pos = m_f.tellp();
        auto sz_of_this_block = current_pos - pos_to_write_in_offset;
        assert(sz_of_this_block > 0);
        //m_f.seekp( pos_to_write_in_offset ); assert(read_uint64_t() == 0xffffffffffffffffuL );
        m_f.seekp( pos_to_write_in_offset ); save_uint64_t( (uint64_t)sz_of_this_block );
        m_f.seekp(current_pos); // move back to here before returning
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

        // New method, specially tuned to GT fields
        zlib_vector:: vec_t doubly_compressed =
            zlib_vector:: deflate(
                vcfGTz:: special_encoder_for_list_of_GT_fields:: deflate(
                    just_last_fields
                )
            );
        save_vector_of_char_with_leading_size( doubly_compressed );
        //auto for_verification = vcfGTz:: special_encoder_for_list_of_GT_fields:: inflate( zlib_vector:: inflate( doubly_compressed ) ); assert(for_verification ==  (just_last_fields |action:: collect) );


    }

    void save_smart_string0(    string const & s) {
        save_code(vcfGTz:: vcfGTz_codes:: code_null_term_string);
        save_string0(s.c_str());
    }

    void            save_string0(char const *s) {
        m_f << s;
        m_f << '\0';
    }
    void            save_vector_of_char_with_leading_size( zlib_vector:: vec_t const & compressed ) {
        save_code(vcfGTz:: vcfGTz_codes:: code_vector_of_char_with_leading_size_plain_GT);
        uint32_t sz = compressed.size();
        save_uint32_t(sz);

        static_assert(std::is_same< unsigned char const * , decltype(compressed.data()) >{} ,"");
        m_f.write   (   reinterpret_cast<char const *>(compressed.data())
                    ,   compressed.size()
                    );
    }


    template<typename T>
    void            save_uint32_t(T sz) {
        static_assert(std:: is_same<T, uint32_t>{}, "");
        auto arr4bytes = bit_conversions:: convert_to_some_bytes<4>(sz);
        assert(sz == bit_conversions:: convert_from_some_bytes(arr4bytes));
        for(uint8_t u8: arr4bytes) {
            m_f << (char) u8;
        }
    }
    void            save_uint64_t(uint64_t sz) {
        auto arr8bytes = bit_conversions:: convert_to_some_bytes<8>(sz);
        assert(sz == bit_conversions:: convert_from_some_bytes(arr8bytes));
        for(uint8_t u8: arr8bytes) {
            m_f << (char) u8;
        }
    }

    void            save_code(vcfGTz:: vcfGTz_codes code) {
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
    (f.rdbuf() && (f.rdbuf()->is_open())) || DIE("Can't find file [" << arg_input_filename << ']');

    auto r = zip(from::ifstream( f ), range:: ints()) ;

    while(!r.empty() && std::get<0>(range:: front_val(r)).substr(0,2) == "##")
        advance(r);

    int64_t SNP_counter = -2; // because -1 is the 'special' SNP for the header

    vcfGTz_writer writer{arg_output_filename};

    auto remember_start_of_this_block = writer.start_a_new_block();
    writer.save_string0("manylines:GTonly:zlib");

    vector<int64_t> offsets_within_block_for_each_line;

    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        auto && s = x |stdget0;

        ++SNP_counter;

        auto fields = tokenize(s, '\t');

        if(SNP_counter % 10000 == 0) { // that's about once a minute on one of the big files
            PP(ELAPSED(), SNP_counter, fields.at(0), fields.at(1), fields.at(2));
        }

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
        offsets_within_block_for_each_line.push_back( writer.m_f.tellp() - remember_start_of_this_block );
        writer.save_this_line(fields);
    }
    writer.close_this_block(remember_start_of_this_block);

    {
        auto remember_start_of_second_block = writer.start_a_new_block();
        writer.save_string0("offsets.into.previous.block");
        writer.save_uint64_t( offsets_within_block_for_each_line.size() );
        for(uint64_t off : offsets_within_block_for_each_line) {
            writer.save_uint64_t( off );
        }

        writer.close_this_block(remember_start_of_second_block);
    }

    // Normally, a block begins with a (non-zero number) that records its size.
    // Hence, I'll just store '0' to record that the blocks are finished
    // (and therefore [presumably] EOF).
    writer.save_uint64_t(0UL);
}
