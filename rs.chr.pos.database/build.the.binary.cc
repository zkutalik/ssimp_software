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
        if(three_builds_for_all_interesting_snps.size() % 1000000 == 0) {
            PP(three_builds_for_all_interesting_snps.size(), utils:: ELAPSED());
        }
        if(three_builds_for_all_interesting_snps.size() > 1000)
            break;
    }

    PP(three_builds_for_all_interesting_snps.size(), utils:: ELAPSED());

    string line;

    ifstream  fs_1kg(argv[1]); // first argument, i.e. should be  "liftover/liftover_hg18.bed"
    int line_number = 0;
    while(getline( fs_1kg, line)) {
        ++line_number;
        if(line_number % 1'000'000 == 0) {
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
        assert(split.at(3).substr(0,2) == "rs");
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
            assert(entry->second.chrom == -1);
            assert(entry->second.hg18  == -1);
            assert(entry->second.hg19  == -1);
            assert(entry->second.hg20  == -1);
            entry->second.chrom = chrom;
            entry->second.hg18  = posa; // zero-based, I think - https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            assert(entry->second.chrom != -1);
            assert(entry->second.hg18  != -1);
            assert(entry->second.hg19  == -1);
            assert(entry->second.hg20  == -1);
        }
    }

    return 0;

    set<int> all_rs_in_any_file;

    ifstream  fs_hrc(argv[2]);
    while(getline( fs_hrc, line)) {
        auto rs = utils:: lexical_cast<int>(line);
        //PP(line, rs);
        all_rs_in_any_file.insert(rs);
    }

    PP(all_rs_in_any_file.size());

    ifstream  fs_uk10k(argv[3]);
    while(getline( fs_uk10k, line)) {
        if(line == "chr,pos.b19,rs,alts,ref,aaf")
            continue;
        try {
            auto split_on_commas = utils:: tokenize(line, ',');
            auto rs_string = split_on_commas.at(2);
            rs_string.substr(0,2) == "rs" || DIE("uk10k, rs...");
            rs_string = rs_string.substr(2);
            auto rs = utils:: lexical_cast<int>(rs_string);
            //PP(line, rs_string, rs);
            all_rs_in_any_file.insert(rs);
        } catch(std::invalid_argument & e) {
            PP(line);
            DIE("uk10k");
        }
    }

    PP(all_rs_in_any_file.size());

    fs_just_rss.eof() || DIE("fs_just_rss.eof()");    // 81'237'339
    fs_1kg.eof() || DIE("fs_1kg.eof()");    // 81'237'339
    fs_hrc.eof() || DIE("fs_hrc.eof()");    // 36'117'264
    fs_uk10k.eof() || DIE("fs_uk10k.eof()");// 17'769'027

    // 88'047'858 in the union

}
