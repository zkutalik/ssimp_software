// Run this on hpc1 to combine everything together:
//
// /data/sgg/sina/public.data/dbsnp
// head dbsnp_hg20.bed liftover_hg18.bed liftover_hg19.bed
//
// chr13   31872704        31872705        rs3     13
// chr13   31873084        31873085        rs4     13
// chr7    92209795        92209796        rs5     7
// chr7    92117816        92117817        rs6     7
// chr7    92150242        92150243        rs7     7

#include<fstream>
#include<cassert>

#include"bits.and.pieces/PP.hh"
#include"bits.and.pieces/utils.hh"
#include"format/format.hh"

using std:: ifstream;
using std:: string;
using utils:: operator <<;

constexpr int CHR_X = 23;
constexpr int CHR_Y = 24;
constexpr int CHR_MT = 25;
constexpr int CHR_PAR = 26; // is this a real thing?

struct reader {
    int     m_build;
    ifstream m_bed;
    int     m_linenumber = 0;

    int m_rs;
    int m_chr;
    int m_pos;

    reader(int build, const char * fn) : m_build(build), m_bed(fn) {
        advance();
    }

    void advance() {
        ++m_linenumber;
        assert(m_rs != -1234);
        string line;
        assert(m_bed);
        getline(m_bed, line);
        if(!m_bed) {
            PP(__LINE__, m_rs,m_chr,m_pos);
            m_rs = -1234;
            assert(m_bed.eof());
            return;
        }

        auto tks = utils:: tokenize(line, '\t');
        assert(tks.size() == 5);

        auto const & rs_str = tks.at(3);
        assert(rs_str.substr(0,2) == "rs");
        int pos0= utils::lexical_cast<int> (tks.at(1));

        m_rs  = utils::lexical_cast<int>(rs_str.substr(2,string::npos));
        assert(tks.at(0).substr(0,3)=="chr");
        auto chr_name = tks.at(0).substr(3,
                            tks.at(0).find('_')-3
                        );
        if(chr_name == "Un") {
            advance();
            return;
        }
        m_chr =     chr_name  == "X" ? CHR_X
                :   chr_name  == "Y" ? CHR_Y
                :   chr_name  == "MT" ? CHR_MT
                :   chr_name  == "PAR" ? CHR_PAR
                :   utils::lexical_cast<int>( chr_name );
        m_pos = utils::lexical_cast<int> (tks.at(2));

        assert(pos0+1 == m_pos);
    }

    bool empty() {
        return m_rs == -1234; // special code to mark the end
    }

};

int main() {
    std:: cout.imbue(std::locale("")); // user-specified locale

    reader hg20 (20, "/data/sgg/sina/public.data/dbsnp/dbsnp_hg20.bed");
    reader hg19 (19, "/data/sgg/sina/public.data/dbsnp/liftover_hg19.bed");
    reader hg18 (18, "/data/sgg/sina/public.data/dbsnp/liftover_hg18.bed");


    //  18 , 324'006'561 ,   678.24 , 1'057'519'607 , 247'199'685
    //  19 , 324'251'184 , 1'353.52 , 1'057'519'607 , 249'240'618
    //  20 , 324'853'746 ,   712.38 , 1'057'519'607 , 248'946'419
    //
    //  18 , 324'006'561 ,   698    , 1'057'519'607 , 247'199'685
    //  19 , 324'251'184 , 1'400.92 , 1'057'519'607 , 249'240'618
    //  20 , 324'853'746 , 1'058.67 , 1'057'519'607 , 248'946'419

    for(auto one_bed_ptr : { &hg20, &hg18, &hg19 }) {
        reader & x = *one_bed_ptr;
        PP(x.m_build, utils:: ELAPSED());
        int max_rs = 0; // TODO: initialize to zero and rerun
        int max_pos = 0; // TODO: initialize to zero and rerun
        while(!x.empty()) {
            assert(x.m_rs > 0);
            assert(x.m_pos > 0);
            max_rs  = std::max(max_rs , x.m_rs );
            max_pos = std::max(max_pos, x.m_pos);
            if(x.m_linenumber % 10000000 == 0) {
                PP(x.m_linenumber, utils:: ELAPSED(), max_rs, max_pos);
            }

            x.advance();
        }
        PP(x.m_build, x.m_linenumber, utils:: ELAPSED(), max_rs, max_pos);
        std:: cout << '\n';

        // 31 bits for the rs number
        //  5 bits for the chr
        // 28 for pos
        // exactly 64 altogether!
    }
}

