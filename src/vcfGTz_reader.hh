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
        std:: cout << "cd = " << (int) cd << '\n';
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
