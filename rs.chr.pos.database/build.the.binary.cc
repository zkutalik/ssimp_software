#include <iostream>
#include <fstream>
#include <set>
#include <map>

#include "../src/bits.and.pieces/PP.hh"
#include "../src/bits.and.pieces/DIE.hh"
#include "../src/bits.and.pieces/utils.hh"

using std:: ifstream;
using std:: cout;
using std:: string;
using std:: set;
using std:: map;
using utils:: operator<<;

struct three_builds_t {
    enum class XYZ : int {
        CHROM_Z = -1,
        CHROM_Y = -2,
        CHROM_X = -3,
    };
    int32_t chrom;
    int32_t hg18;
    int32_t hg19;
    int32_t hg20;

    three_builds_t() : chrom(-1), hg18(-1), hg19(-1), hg20(-1) {}
};

int main(int argc, char **argv) {
    argc == 6 || DIE("needs 6 args");

    map<int, three_builds_t> three_builds_for_all_interesting_snps;

    ifstream fs_just_rss(argv[4]);
    int32_t one_rs_number;
    while(fs_just_rss >> one_rs_number) {
        assert(three_builds_for_all_interesting_snps.count(one_rs_number)==0);
        three_builds_for_all_interesting_snps[one_rs_number];
        assert(three_builds_for_all_interesting_snps.count(one_rs_number)==1);
        if(three_builds_for_all_interesting_snps.size() % 10'000'000 == 0) {
            PP(three_builds_for_all_interesting_snps.size(), utils:: ELAPSED());
        }
        //if(three_builds_for_all_interesting_snps.size() > 1000) break;
    }

    PP(three_builds_for_all_interesting_snps.size(), utils:: ELAPSED());
    assert(three_builds_for_all_interesting_snps.size() == 88'047'858);

    for(int argv_index : {1,2,3}) { // should be filenames to hg18, hg19, and hg20
        PP(argv_index, argv[argv_index], utils:: ELAPSED());
        ifstream  fs_one_build(argv[argv_index]);

        int line_number = 0;
        string line;

        while(getline( fs_one_build, line)) {
            ++line_number;
            if(line_number % 10'000'000 == 0) {
                PP(line_number, utils:: ELAPSED());
            }
            auto split = utils:: tokenize(line, '\t');
            //PP(line, split);
            auto chrom_string = split.at(0);
            assert(chrom_string.substr(0,3) == "chr");
            chrom_string = chrom_string.substr(3);
            if(chrom_string.size() >= 7) {
                if(chrom_string.substr( chrom_string.size() - 7 ) == "_random")
                    continue; // Skip the _random ones. TODO: Confirm this is OK
            }
            if(chrom_string.substr(0,3) == "Un_")
                continue;
            if(chrom_string             == "PAR")
                continue;
            if(chrom_string             == "MT")
                continue;
            assert(split.at(3).substr(0,2) == "rs");

            try {
                int chrom = [&]() ->int {
                    if(chrom_string=="X") return static_cast<int>( three_builds_t:: XYZ:: CHROM_X );
                    if(chrom_string=="Y") return static_cast<int>( three_builds_t:: XYZ:: CHROM_Y );
                    return utils:: lexical_cast<int>(chrom_string);
                }();
                int posa  = utils:: lexical_cast<int>(split.at(1));
                int posb  = utils:: lexical_cast<int>(split.at(2));
                int rs    = utils:: lexical_cast<int>(split.at(3).substr(2));
                //PP(chrom, posa, rs);
                assert(posb == posa+1);

                auto entry = three_builds_for_all_interesting_snps.find(rs);
                if(entry != three_builds_for_all_interesting_snps.end()) {
                    if(entry->second.chrom == -1) {
                        entry->second.chrom = chrom;
                    }
                    if(entry->second.chrom != chrom) {
                        PP(line, split, chrom, entry->second.chrom);
                    }
                    assert(entry->second.chrom == chrom);

                    int * position_under_current_build = nullptr;
                    switch(argv_index) {
                        break; case 1: position_under_current_build = & entry->second.hg18;
                        break; case 2: position_under_current_build = & entry->second.hg19;
                        break; case 3: position_under_current_build = & entry->second.hg20;
                        break; default: DIE("argv_index out of range");
                    }
                    assert(position_under_current_build);
                    assert(*position_under_current_build  == -1);
                    *position_under_current_build  = posa; // zero-based, I think - https://genome.ucsc.edu/FAQ/FAQformat.html#format1
                    assert(*position_under_current_build  != -1);
                }
            } catch (std:: invalid_argument &) {
                PP(line, split);
                exit(1);
            }
        }
        fs_one_build.eof() || DIE("fs_one_build.eof()");
    }

}
