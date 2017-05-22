
#include "libStatGen/include/VcfFileReader.h"

#include "bits.and.pieces/utils.hh"
#include "bits.and.pieces/PP.hh"

void first_attempt_at_vcfgztbi_file() {
    VcfFileReader reader;
    VcfHeader header;

    // reader.open will throw if the file doesn't exist
    auto ret = reader.open("/data/sgg/aaron/shared/ref_panels/1kg/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", header);

    auto ret2 = reader.readVcfIndex(); // should find the .tbi file nearby

    utils:: print_type(ret);
    utils:: print_type(ret2);
    PP(ret);
    PP(ret2);
    assert(ret);
    assert(ret2);
    reader.set1BasedReadSection("22", 18'000'000, 19'100'000);
    reader.setDiscardRules(     VcfFileReader:: DISCARD_FILTERED
                            |   VcfFileReader:: DISCARD_NON_PHASED
                            |   VcfFileReader:: DISCARD_MISSING_GT
                            |   VcfFileReader:: DISCARD_MULTIPLE_ALTS
                            );
    reader.addDiscardMinMinorAlleleCount(100000, NULL);

    VcfRecord record;
    while(reader.readRecord(record))
    {
        // Your record specific processing here.
        if(1)
        PP( record.getChromStr()
          , record.get1BasedPosition()
          , record.getIDStr()
          , record.getNumSamples()         // TODO: Check this is constant?
          , record.hasAllGenotypeAlleles() // TODO: Should always be true?
          );
        //PP( record.getNumGTs(0) , record.getNumGTs(1) , record.getNumGTs(2) , record.getGT(0, 0) , record.getGT(0, 1) , record.getGT(6, 0) , record.getGT(6, 1));


    }
        PP(reader.getNumRecords());
        PP(reader.getNumKeptRecords());
        PP(reader.getTotalReadRecords());
    reader.close();
    exit(0);
}
