#include <string>

//TODO: move more of this into the .cc file, especially the #includes
#include "libStatGen/include/VcfFileReader.h"
#include "format/format.hh"
#include "file.reading.hh" // just for 'chrpos', I think
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"

namespace tbi {
    struct read_vcf_with_tbi {
        VcfFileReader reader;
        VcfHeader header;

        read_vcf_with_tbi(std:: string filename) {
            auto ret = reader.open( filename.c_str() , header);

            auto ret2 = reader.readVcfIndex(); // should find the .tbi file nearby

            ret || DIE("Can't open that vcf(.gz) file: [" << filename << "]");
            ret2 || DIE("Can't open that vcf(.gz).tbi file: [" << filename << "]");

            reader.setDiscardRules  (   VcfFileReader:: DISCARD_FILTERED
                                    |   VcfFileReader:: DISCARD_NON_PHASED
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
    };
}
void first_attempt_at_vcfgztbi_file();
