#include<iostream>

#include "gzstream.h"

#include "range/range_from.hh"
#include "range/range_view.hh"
#include "range/range_action.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "zlib-vector-of-char/zlib-vector.hh"

namespace from = range:: from;
namespace action = range:: action;

using std:: vector;
using utils:: tokenize;
using utils:: stdget0;
using utils:: stdget1;
using range:: range_from_begin_end;
using range:: view:: take;

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

    int64_t SNP_counter = -2; // because -1 is the 'special' SNP for the header

    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        (void)line_no;
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
        auto INFO = fields.at(7);
        auto FORMAT = fields.at(8);
        fields.erase( fields.begin() + 7, fields.begin() + 9 );

        { // This is pointless, just an assertion to help check my compression
            zlib_vector:: vec_t as_a_vector{ s.begin(), s.end() };
            assert(s.size() == as_a_vector.size());
            auto compressed = zlib_vector:: deflate( as_a_vector );
            auto uncompressed = zlib_vector:: inflate( compressed );
            PP(as_a_vector.size(), compressed.size());
            assert(uncompressed == as_a_vector);
        }
    }
}
