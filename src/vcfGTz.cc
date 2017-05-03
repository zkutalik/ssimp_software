#include<iostream>

#include "gzstream.h"

#include "range/range_from.hh"
#include "range/range_action.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "zlib-vector-of-char/zlib-vector.hh"

namespace from = range:: from;
namespace action = range:: action;

using std:: vector;
using utils:: stdget0;
using utils:: stdget1;

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

    gz:: igzstream f(argv[1]);
    f.fail() && DIE("Can't find file [" << argv[1] << ']');

    auto r = zip(from::ifstream( f ), range:: ints()) ;

    while(!r.empty() && std::get<0>(range:: front_val(r)).substr(0,2) == "##")
        advance(r);

    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        (void)line_no;
        auto && s = x |stdget0;

        //PP(line_no, s);
        zlib_vector:: vec_t as_a_vector{ s.begin(), s.end() };
        assert(s.size() == as_a_vector.size());
        auto compressed = zlib_vector:: deflate( as_a_vector );
        { // This is pointless, just an assertion to help check my compression
            auto uncompressed = zlib_vector:: inflate( compressed );
            PP(as_a_vector.size(), compressed.size());
            assert(uncompressed == as_a_vector);
        }
    }
}
