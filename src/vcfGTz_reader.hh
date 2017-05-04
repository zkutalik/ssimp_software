#include "bits.and.pieces/utils.hh"

#include<fstream>
#include<array>

using std:: ios_base;
using utils:: ssize;

namespace vcfGTz {
    enum class vcfGTz_codes : uint8_t {
                code_end_of_file        = 0
            ,   code_null_term_string   = 1 // null-terminated string - when everything else fails
            ,   code_vector_of_char_with_leading_size_plain_GT = 2 // for storing the zlib-compressed stuff
    };

    namespace bit_conversions {
    template<size_t N, typename T>
    static
    std:: array<uint8_t,N> convert_to_some_bytes(T u) {
        static_assert(N==8 || N==4, "");
        static_assert(      ( std:: is_same<T, uint32_t>{} && N==4 )
                         || ( std:: is_same<T, uint64_t>{} && N==8 ), "");

        std:: array<uint8_t,N> arr;
        static_assert(N==utils::ssize(arr) ,"");
        for(int byte = 0; byte < utils::ssize(arr); ++byte) {
            auto by = (u >> (byte*8)) & 255;
            arr.at(N-1-byte) = by;
        }
        return arr;
    }
    template<size_t N>
    static
    auto convert_from_some_bytes(std:: array<uint8_t,N> arr) {
        static_assert(N == 4 || N==8 ,""); // to return a uint32_t
        using return_type = std::conditional_t<N==4, uint32_t, uint64_t>;
        static_assert(sizeof(return_type) == N ,"");
        return_type u = 0;
        for(int byte = arr.size(); byte > 0; --byte) {
            u <<= 8;
            u += arr.at(N-byte);
        }
        return u;
    }
    }; // namespace bit_conversions

struct special_encoder_for_list_of_GT_fields {
    /* for more efficiency, instead of simply concatenating
     * the fields with tabs and then using zlib, I use
     * this more customized method.
     *
     * Note, this is lossless, i.e, it works in the header line
     * even though it doesn't have entries like '0|0'
     *
     * HOWEVER, this assumes (reasonably!) that the original
     * strings don't have '\0' in them. This is what I
     * call ensure_there_are_no_nulls() before processing
     * any line.
     */
    using vuc_t = std:: vector<unsigned char>; // same as 'zlib_vector:: vec_t

    template<typename T> // T will be a range-type
    static
    vuc_t   deflate(T just_last_fields) {
        vuc_t special_encoding_of_list_of_GT_fields;
        for(auto & one_GT_field : just_last_fields) {
            if  (   one_GT_field.size() == 3
                 && (one_GT_field.at(1) == '|' || one_GT_field.at(1) == '/' )
                 && (one_GT_field.at(0) >= '0' || one_GT_field.at(0) <= '3' )
                 && (one_GT_field.at(2) >= '0' || one_GT_field.at(2) <= '3' )
                ) {
                static_assert('@' == 64 ,"not ASCII?");
                // encode this in one character
                bool slash_not_pipe = one_GT_field.at(1) == '/';
                int left  = one_GT_field.at(0) - '0';
                int right = one_GT_field.at(2) - '0';
                assert(left  >= 0 && left  <= 3); // i.e. just two bits
                assert(right >= 0 && right <= 3); // i.e. just two bits
                char together = '@'
                                | slash_not_pipe
                                | (left << 1)
                                | (right << 3)
                                ;
                assert(together >= '@' && together <= '_');
                special_encoding_of_list_of_GT_fields.push_back( together );
            }
            else {
                special_encoding_of_list_of_GT_fields.push_back( '\t' );
                special_encoding_of_list_of_GT_fields.insert(   special_encoding_of_list_of_GT_fields.end()
                                                            ,   one_GT_field.begin()
                                                            ,   one_GT_field.end()
                                                            );
                special_encoding_of_list_of_GT_fields.push_back( '\0' );
            }
        }
        return special_encoding_of_list_of_GT_fields;
    }

    static
    std::vector<std::string>   inflate(vuc_t const & special_encoding_of_list_of_GT_fields) {
        std::vector<std:: string> vector_of_GT_fields;
        auto it = special_encoding_of_list_of_GT_fields.begin();
        auto end = special_encoding_of_list_of_GT_fields.end();
        while(it != end) {
            std:: string one_GT_field;

            char typ = *it;
            ++it;

            if(typ == '\t') {
                assert(it != end);
                while(*it != '\0') {
                    one_GT_field.push_back(*it);
                    ++it;
                    assert(it != end);
                }
                assert(it != end);
                assert(*it == '\0');
                ++it;
            }
            else {
                assert( typ >= '@' && typ <= '_');
                int together = typ - '@';
                bool    slash_not_pipe  = together & 1;
                int     left            = (together>>1) & 3;
                int     right           = (together>>3) & 3;
                one_GT_field.push_back('0' + left );
                one_GT_field.push_back(slash_not_pipe ? '/' : '|');
                one_GT_field.push_back('0' + right);
            }
            vector_of_GT_fields.push_back(one_GT_field);
        }
        return vector_of_GT_fields;
    }
};

struct vcfGTz_reader {
    std:: ifstream m_f;
    decltype(m_f.tellg())   m_my_beginning;
    vcfGTz_reader(std:: string file_name) : m_f {   file_name
                                                ,   ios_base::in
                                                  | ios_base::binary
    }{
        m_my_beginning = m_f.tellg();
    }

    std:: string read_smart_string0() {
        auto cd = read_code();
        assert( cd == vcfGTz:: vcfGTz_codes:: code_null_term_string );
        return read_string0();
    }
    std:: string read_string0() {
        std:: string s;
        getline(m_f, s, '\0');
        return s;
    }
    auto            read_vector_of_char_with_leading_size() {
        auto cd = read_code();
        assert(cd == vcfGTz_codes:: code_vector_of_char_with_leading_size_plain_GT);
        auto length = read_uint32_t();
        std:: vector<unsigned char> vec(length);
        m_f.read( reinterpret_cast<char*>(vec.data()), length );
        return vec;
    }
    vcfGTz_codes    read_code() {
        assert(m_f);
        uint8_t u = (char) m_f.get();
        assert(m_f);
        return static_cast<vcfGTz_codes>(u);
    }

    int64_t     read_offset_at_start_of_block();

    uint32_t    read_uint32_t();
    uint64_t    read_uint64_t();
};

} // namespace vcfGTz
