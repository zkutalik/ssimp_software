#include "vcfGTz_reader.hh"

#include "bits.and.pieces/utils.hh"
#include "bits.and.pieces/PP.hh"
#include "range/range.hh"

#include<cassert>

using utils:: ssize;
using utils:: save_ostream_briefly;
using std:: cout;

namespace vcfGTz {
    int64_t     vcfGTz_reader:: read_offset_at_start_of_block() {
        int64_t i = read_uint64_t();
        PP(i);
        assert(i>=0);
        assert(i> 0);
        return i;
    }
    uint64_t    vcfGTz_reader:: read_uint64_t() {
        constexpr int N = 8;
        std:: array<uint8_t, N> arr;
        static_assert(N == sizeof(this->read_uint64_t()) ,"");
        for(int i=0; i<ssize(arr); ++i) {
            assert(m_f);
            m_f >> arr.at(i);
        }
        assert(m_f);
        return   bit_conversions:: convert_from_some_bytes(arr);
    }
} // namespace vcfGTz
