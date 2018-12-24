#include"file.reading.vcfgztbi.hh"

#include "bits.and.pieces/utils.hh"
#include "bits.and.pieces/PP.hh"
#include "bits.and.pieces/DIE.hh"
#include "format/format.hh"

#include "options.hh"

#include <numeric>
#include <fstream>

constexpr double THRESHOLD_OF_ACCEPTABLE_REF_MISSINGNESS = 0.05;

namespace tbi {
        read_vcf_with_tbi :: read_vcf_with_tbi(std:: string filename, int chromosome) {
            // First, test if 'filename' exists. If not, replace any embedded '{CHRM}' with the chromosome number
            {
                std:: ifstream test_if_file_exists(filename.c_str());
                if(!test_if_file_exists) {
                    auto CHRM = filename.find("{CHRM}");
                    if(chromosome == 23)
                    {
                        if(filename.substr(CHRM) == "{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
                        {
                            filename = filename.substr(0, CHRM) + "{CHRM}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz";
                        }
                    }
                    std:: ostringstream filename_for_this_chromosome;
                    filename_for_this_chromosome << filename.substr(0, CHRM);
                    if(chromosome == 23)
                        filename_for_this_chromosome << "X";
                    else
                        filename_for_this_chromosome << chromosome;
                    filename_for_this_chromosome << filename.substr(CHRM+6); // 6, to skip over {CHRM}
                    filename = filename_for_this_chromosome.str();
                }
            }
            auto ret = [&](){
                if(options:: opt_sample_names.empty()) {
                    return reader.open( filename.c_str() , header );
                }
                else {
                    return reader.open( filename.c_str() , header, options:: opt_sample_names.c_str(), NULL, NULL );
                }
            }();

            auto ret2 = reader.readVcfIndex(); // should find the .tbi file nearby

            ret || DIE("Can't open that vcf(.gz) file: [" << filename << "]");
            ret2 || DIE("Can't open that vcf(.gz).tbi file: [" << filename << "]");

            reader.setDiscardRules  (   0
                                    //|   VcfFileReader:: DISCARD_FILTERED // don't use this after all.
                                    //|   VcfFileReader:: DISCARD_NON_PHASED // TODO: Ask SR and ZK if we should include this. This discards '/' doesn't it?
                                    //|   VcfFileReader:: DISCARD_MISSING_GT // As of December 2018, we don't discard SNVs with some missing samples. Instead we replace them with the mean
                                    |   VcfFileReader:: DISCARD_MULTIPLE_ALTS
                                    );
            reader.addDiscardMinMinorAlleleCount(1, NULL); // otherwise, there is no variation and it's not very useful
        }
        void read_vcf_with_tbi:: set_region(file_reading:: chrpos b, file_reading:: chrpos e) {
            assert(b.chr == e.chr);
            if(b.chr == 23)
                reader.set1BasedReadSection("X"                                       , b.pos, e.pos);
            else
                reader.set1BasedReadSection(AMD_FORMATTED_STRING("{0}", b.chr).c_str(), b.pos, e.pos);
            /* "...will set the code to read the specified chromosome starting
             * at the specified 1-based position up to, but not including, the
             * specified 1-based end position"
             * This is what we want - http://genome.sph.umich.edu/wiki/LibStatGen:_VCF#Specifying_Discard_Rules
             */
        }
    double af_under_missingness(VcfRecord & record) {
        // Computes the allele frequency for this reference SNV.
        //
        // the average of the zeroes and ones, i.e. ignoring the cells that
        // contain -2 (missing data) or -1 (missing X chromosome for male)
        int const N = record.getNumSamples();
        double total_g = 0;
        int count_g = 0;
        for(int i=0;i<N;++i) {
            int l = record.getGT(i, 0);
            int r = record.getGT(i, 1);
            if(l>=0) {
                total_g += l;
                count_g += 1;
            }
            if(r>=0) {
                total_g += r;
                count_g += 1;
            }
        }
        return total_g / count_g;
    }
    RefRecord   convert_VcfRecord_to_RefRecord(VcfRecord & record) {
        RefRecord rr;
        rr.pos      =   record.get1BasedPosition();
        rr.ID       =   record.getIDStr();
        rr.ref      =   record.getRefStr();
        rr.alt      =   record.getAltStr();
        assert(record.getNumAlts() == 1); // The 'DISCARD' rule should already have skipped those with more alts
        auto maf = af_under_missingness(record);

        int const N = record.getNumSamples(); // TODO: verify this is the same in every SNP?
        int missing_in_the_reference_panel = 0;
        rr.z12.reserve(N);
        for(int i=0;i<N;++i) {
            //assert(2==record.getNumGTs(i));
            int l = record.getGT(i, 0);
            int r = record.getGT(i, 1);
            if(r==-1)
            {
                // This is a male on the X chromosome.
                // We'll upscale by copying the left into the right.
                // i.e.  0/-1 becomes 0/0
                //  and  1/-1 becomes 1/1
                r = l;
            }
            if(l==-2) // missing in the reference panel - e.g. ".|." in the vcf file
                ++missing_in_the_reference_panel;
            if(l==-2) // missing in the reference panel - e.g. ".|." in the vcf file
                l = maf;
            if(r==-2) // missing in the reference panel - e.g. ".|." in the vcf file
                r = maf;
            rr.z12.push_back(l+r);
        }
        assert(N == utils:: ssize(rr.z12));

        assert(maf >= 0.0);
        assert(maf >  0.0);
        assert(maf <= 1.0);
        assert(maf <  1.0);
        if(maf > 0.5) maf = 1.0-maf;
        assert(maf <= 0.5);
        //assert(maf <  0.5);
        rr.maf = maf;
        rr.proportion_of_missing_ref_data = missing_in_the_reference_panel / double(N);
        return rr;
    }
    bool    read_vcf_with_tbi:: read_record_into_a_RefRecord(RefRecord &rr) {
                /* Note that this function calls itself recursively in some
                 * cases in order to skip a record that isn't desired.
                 * Therefore, after considering the recursion, its semantics
                 * are to either:
                 *  - set rr to a valid record and return true
                 *  OR
                 *  - return false if no such record is found
                 */
                VcfRecord record;
                bool b = reader.readRecord(record);
                if(b && record.getNumAlts() != 1) {
                    // On some machines, DISCARD_MULTIPLE_ALTS doesn't work,
                    // so I'll manually do it here
                    return read_record_into_a_RefRecord(rr);
                }
                if(b) {
                    rr = tbi:: convert_VcfRecord_to_RefRecord(record);
                    if(rr.proportion_of_missing_ref_data > THRESHOLD_OF_ACCEPTABLE_REF_MISSINGNESS)
                        return read_record_into_a_RefRecord(rr);
                }
                return b;
    }
} // namespace tbi
