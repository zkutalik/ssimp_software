#include <string>

//TODO: move more of this into the .cc file, especially the #includes
#include "libStatGen/include/VcfFileReader.h"
#include "file.reading.hh" // just for 'chrpos', I think

#include<algorithm>

namespace tbi {
    struct RefRecord;
    struct read_vcf_with_tbi {
        VcfFileReader reader;
        VcfHeader header;

        read_vcf_with_tbi(std:: string filename, int chromosome); // 'chromosome' to take the place of '{CHRM}' in the filename

        void set_region(file_reading:: chrpos b, file_reading:: chrpos e);
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
    RefRecord   convert_VcfRecord_to_RefRecord(VcfRecord & record);
}
