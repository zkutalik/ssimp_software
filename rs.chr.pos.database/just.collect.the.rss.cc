#include <iostream>
#include <fstream>
#include <set>

#include "../src/bits.and.pieces/PP.hh"
#include "../src/bits.and.pieces/DIE.hh"
#include "../src/bits.and.pieces/utils.hh"

using std:: ifstream;
using std:: cout;
using std:: string;
using std:: set;

int main(int argc, char **argv) {
    argc == 5 || DIE("needs 5 args");

    string line;
    set<int> all_rs_in_any_file;

    ifstream  fs_1kg(argv[1]);
    while(getline( fs_1kg, line)) {
        if(line == "#CHROM\tPOS\tID")
            continue;
        //PP(line);
        auto tokens = utils:: tokenize(line, '\t');
        auto third = tokens.at(2);

        for(auto split_on_semicolons : utils:: tokenize(third, ';')) {
            if(split_on_semicolons.substr(0, 2) == "rs") {
                split_on_semicolons = split_on_semicolons.substr(2);
                auto rs_1kg = utils:: lexical_cast<int>(split_on_semicolons);
                //PP(split_on_semicolons); PP(rs_1kg);
                all_rs_in_any_file.insert(rs_1kg);
            }
        }

    }

    PP(all_rs_in_any_file.size());

    ifstream  fs_hrc(argv[2]);
    while(getline( fs_hrc, line)) {
        auto rs = utils:: lexical_cast<int>(line);
        //PP(line, rs);
        all_rs_in_any_file.insert(rs);
    }

    PP(all_rs_in_any_file.size());

    ifstream  fs_uk10k(argv[3]);
    while(getline( fs_uk10k, line)) {
        if(line, line.substr(0,27)  == "chr,pos.b19,rs,alts,ref,aaf")
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

    fs_1kg.eof() || DIE("fs_1kg.eof()");    // 81'237'339
    fs_hrc.eof() || DIE("fs_hrc.eof()");    // 36'117'264
    fs_uk10k.eof() || DIE("fs_uk10k.eof()");// 17'769'027

    // 88'047'858 in the union

    std:: ofstream out(argv[4]);
    for(auto && rs : all_rs_in_any_file) {
        out << rs << '\n';
    }
    out || DIE("problem with output");
    out.close();
    out || DIE("problem with output");
}
