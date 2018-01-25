#include <iostream>
#include <fstream>
#include <set>
#include <map>

#include "../src/bits.and.pieces/PP.hh"
#include "../src/bits.and.pieces/DIE.hh"
#include "../src/bits.and.pieces/utils.hh"

using std:: ifstream;
using std:: ofstream;
using std:: cout;
using std:: string;
using std:: set;
using std:: map;
using utils:: operator<<;

struct three_builds_t {
    enum class XYZ : int {
        CHROM_unknown = -1,
        CHROM_Y = -2,
        CHROM_X = -3,
    };
    int32_t chrom;
    int32_t hg18;
    int32_t hg19;
    int32_t hg20;

    three_builds_t() : chrom(-1), hg18(-1), hg19(-1), hg20(-1) {}
};

void write_int(int32_t & input, ofstream & o) {
    for(int i=3; i>=0; --i) {
        unsigned char one_byte = (input >> (i*8)) & 0xff;
        o.write( reinterpret_cast<const char *>(& one_byte), 1 );
    }
    o || DIE("couldn't write");
}

int main(int argc, char **argv) {
    argc == 6 || DIE("needs 6 args");

    map<int, three_builds_t> three_builds_for_all_interesting_snps;

    ifstream fs_just_rss(argv[4]); // this is the list of approx. 88 million rs-numbers
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
    assert(three_builds_for_all_interesting_snps.size() == 88'047'858
        || three_builds_for_all_interesting_snps.size() == 89'506'168);

    for(int argv_index : {1,2,3}) { // the filenames to hg18, hg19, and hg20
        PP(argv_index, argv[argv_index], utils:: ELAPSED());
        ifstream  fs_one_build(argv[argv_index]);

        int line_number = 0;
        string line;

        while(getline( fs_one_build, line)) { // read a line from the liftover file
            ++line_number;
            if(line_number % 10'000'000 == 0) {
                PP(line_number, utils:: ELAPSED());
            }
            // next few lines parse the relevant info from the liftover line
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

            // now we try to store the data that we have just parsed into
            // the datastructure

            try {
                int chrom = [&]() ->int {
                    if(chrom_string=="X")
                        return static_cast<int>( three_builds_t:: XYZ:: CHROM_X );
                    if(chrom_string=="Y")
                        return static_cast<int>( three_builds_t:: XYZ:: CHROM_Y );
                    int chrom = utils:: lexical_cast<int>(chrom_string);
                    assert(chrom >= 1 && chrom <= 22); // unless it's "X" or "Y", of course.
                    return chrom;
                }();
                int posa  = utils:: lexical_cast<int>(split.at(1));
                int posb  = utils:: lexical_cast<int>(split.at(2));
                int rs    = utils:: lexical_cast<int>(split.at(3).substr(2));
                //PP(chrom, posa, rs);
                assert(posb == posa+1);

                // try to find this snp in the database. The database started
                // with nearly 90 million SNPs, but some might have been discarded
                // already (see below) if two builds assigned the SNP to different chromosomes
                auto entry = three_builds_for_all_interesting_snps.find(rs);
                if(entry != three_builds_for_all_interesting_snps.end()) { // ... if it's in the database

                    // We need to update the entry in the database with the data we've just
                    // read from the liftover file

                    // First, record the chromosome number if not already recorded
                    if(entry->second.chrom == -1) {
                        entry->second.chrom = chrom;
                    }

                    // Second, if the chromosome contradicts that from an earlier build,
                    // then erase this entry from the database and no longer consider this rs-number.
                    if(entry->second.chrom != chrom) {
                        PP(three_builds_for_all_interesting_snps.size(), line, split, chrom, entry->second.chrom);
                        assert(chrom >= 1 || chrom == -3 || chrom == -2);
                        assert(chrom >= 1 || chrom == -3); // Actually, we're not considering Y for now
                        switch(entry->second.chrom) {
                            break; case int(three_builds_t:: XYZ:: CHROM_X): {}
                            break; case int(three_builds_t:: XYZ:: CHROM_Y): {}
                            break; default: {
                                assert(entry->second.chrom >= 1);
                            }
                        }
                        three_builds_for_all_interesting_snps.erase(entry); // discard this rs-number. Too wierd that the chromosomes aren't the same across builds
                        continue;
                    }
                    assert(entry->second.chrom == chrom);

                    // Finally, store the position in the appropriate place
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

    ofstream fs_out(argv[5], std:: ios_base:: binary);
    map<int, int> mask_of_three_builds_frequency;
    map<int, int> chrom_counts;
    for(auto && rs_builds : three_builds_for_all_interesting_snps) {
        auto rs = rs_builds.first;
        auto three_builds = rs_builds.second;
        if(three_builds.chrom == -3)
            three_builds.chrom = 23; // just rename X/-3 to 23 when actually writing the final database

        write_int(rs, fs_out);
        write_int(three_builds.chrom, fs_out);
        write_int(three_builds.hg18, fs_out);
        write_int(three_builds.hg19, fs_out);
        write_int(three_builds.hg20, fs_out);

        auto mask_of_three_builds = (three_builds.hg18 != -1 ? 1 : 0)
                                  + (three_builds.hg19 != -1 ? 2 : 0)
                                  + (three_builds.hg20 != -1 ? 4 : 0)
            ;
        ++mask_of_three_builds_frequency[mask_of_three_builds];
        ++chrom_counts[three_builds.chrom];
    }
    for(auto && chrom_count : chrom_counts)
        PP(chrom_count.first, chrom_count.second);
    fs_out.close();;
    fs_out || DIE("Couldn't close output file");

    for(auto && mask : mask_of_three_builds_frequency) {
        PP(mask.first, mask.second);
    }

}
