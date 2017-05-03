#include<iostream>

#include "gzstream.h"

#include "range/range_from.hh"
#include "range/range_action.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"

namespace from = range:: from;
namespace action = range:: action;

using utils:: stdget0;
using utils:: stdget1;

int main(int argc, char **argv) {
    argc == 3 || DIE(argv[0] << ": requires 2 args, input (.vcf/.vcfgz) and output");

    gz:: igzstream f(argv[1]);
    (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << argv[1] << ']');

    auto r = zip(from::ifstream( f ), range:: ints()) ;
    using utils:: operator<<;
    for(auto && x : r ) {
        int line_no = x |stdget1;
        auto && s = x |stdget0;
        if(line_no == 10)
            break;
        PP(line_no, s);
    }
}
