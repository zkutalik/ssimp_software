#include <string>

//TODO: move more of this into the .cc file, especially the #includes
#include "libStatGen/include/VcfFileReader.h"
#include "format/format.hh"
#include "file.reading.hh" // just for 'chrpos', I think
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"

#include<algorithm>

namespace tbi {
    struct RefRecord;
    struct read_vcf_with_tbi {
        VcfFileReader reader;
        VcfHeader header;

        read_vcf_with_tbi(std:: string filename) {
            auto ret = reader.open( filename.c_str() , header);

            auto ret2 = reader.readVcfIndex(); // should find the .tbi file nearby

            ret || DIE("Can't open that vcf(.gz) file: [" << filename << "]");
            ret2 || DIE("Can't open that vcf(.gz).tbi file: [" << filename << "]");

            reader.setDiscardRules  (   0
                                    |   VcfFileReader:: DISCARD_FILTERED
                                    //|   VcfFileReader:: DISCARD_NON_PHASED // TODO: Ask SR and ZK if we should include this. This discards '/' doesn't it?
                                    |   VcfFileReader:: DISCARD_MISSING_GT
                                    |   VcfFileReader:: DISCARD_MULTIPLE_ALTS
                                    );
            reader.addDiscardMinMinorAlleleCount(1, NULL); // otherwise, there is no variation and it's not very useful
        }

        void set_region(file_reading:: chrpos b, file_reading:: chrpos e) {
            assert(b.chr == e.chr);
            reader.set1BasedReadSection(AMD_FORMATTED_STRING("{0}", b.chr).c_str(), b.pos, e.pos);
            /* "...will set the code to read the specified chromosome starting
             * at the specified 1-based position up to, but not including, the
             * specified 1-based end position"
             * This is what we want - http://genome.sph.umich.edu/wiki/LibStatGen:_VCF#Specifying_Discard_Rules
             */
        }
        inline
        bool    read_record_into_a_RefRecord(RefRecord &rr);
    };
    struct RefRecord {
        int     pos;
        std::string  ID;
        std::string  ref;
        std::string  alt;
        std::vector<int> z12;
        double          maf; // *Minor allele frequency, i.e. always <= 0.50

        bool operator< (file_reading:: chrpos crps) const {
            return pos < crps.pos;
        }
    };
    inline
    bool operator<( file_reading:: chrpos const & crps, RefRecord const & rr) {
        return crps.pos < rr.pos;
    }
    inline
    RefRecord   convert_VcfRecord_to_RefRecord(VcfRecord & record) {
        RefRecord rr;
        rr.pos      =   record.get1BasedPosition();
        rr.ID       =   record.getIDStr();
        rr.ref      =   record.getRefStr();
        rr.alt      =   record.getAltStr();
        assert(record.getNumAlts() == 1); // The 'DISCARD' rule should already have skipped those with more alts
        assert(record.hasAllGenotypeAlleles());
        int const N = record.getNumSamples(); // TODO: verify this is the same in every SNP?
        for(int i=0;i<N;++i) {
            assert(2==record.getNumGTs(i));
            int l = record.getGT(i, 0);
            int r = record.getGT(i, 1);
            assert((l | 1) == 1);
            assert((r | 1) == 1);
            rr.z12.push_back(l+r);
        }
        assert(N == utils:: ssize(rr.z12));
        int total_alternative_alleles = std:: accumulate( rr.z12.begin(), rr.z12.end(), 0);
        double maf = double(total_alternative_alleles) / (2*N);
        //PP(rr.ID, maf);
        assert(maf >= 0.0);
        assert(maf >  0.0);
        assert(maf <= 1.0);
        assert(maf <  1.0);
        if(maf > 0.5) maf = 1.0-maf;
        assert(maf <= 0.5);
        //assert(maf <  0.5);
        rr.maf = maf;
        return rr;
    }
    inline
    bool    read_vcf_with_tbi:: read_record_into_a_RefRecord(RefRecord &rr) {
                VcfRecord record;
                bool b = reader.readRecord(record);
                if(b) {
                    rr = tbi:: convert_VcfRecord_to_RefRecord(record);
                }
                return b;
    }
}
