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
    argc == 3 || DIE(argv[0] << ": requires 2 args, input (.vcf/.vcfgz) and output");

    gz:: igzstream f(argv[1]);
    (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << argv[1] << ']');

    auto r = zip(from::ifstream( f ), range:: ints()) ;

    while(!r.empty() && std::get<0>(range:: front_val(r)).substr(0,2) == "##")
        advance(r);

    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        auto && s = x |stdget0;

        PP(line_no, s);
        zlib_vector:: vec_t as_a_vector{ s.begin(), s.end() };
        assert(s.size() == as_a_vector.size());
        auto compressed = zlib_vector:: deflate( std::move(as_a_vector) );
        {
            auto uncompressed = zlib_vector:: inflate( compressed );
            PP(compressed.size(), uncompressed.size());
        }
    }
}
