#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

#include <gzstream.h>

#include "options.hh"
#include "file.reading.hh"

#include "bits.and.pieces/utils.hh"
#include "mvn/mvn.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "range/range_view.hh"
#include "range/range_from.hh"
#include "range/range_action.hh"
#include "format/format.hh"

#include "file.reading.vcfgztbi.hh"

namespace view = range:: view;
namespace from = range:: from;
namespace action = range:: action;

using std:: cout;
using std:: endl;
using std:: string;
using std:: vector;
using std:: setw;
using std:: unique_ptr;
using std:: make_unique;
using std:: ofstream;
using std:: ostream;

using file_reading:: chrpos;
using file_reading:: SNPiterator;
using file_reading:: GenotypeFileHandle;
using file_reading:: GwasFileHandle;

using utils:: ssize;
using utils:: print_type;
using utils:: stdget0;
using utils:: stdget1;


namespace ssimp {
// Some forward declarations
static
void impute_all_the_regions( file_reading:: GenotypeFileHandle         raw_ref_file
                             , file_reading:: GwasFileHandle             gwas
                             );
static
std:: unordered_map<string, chrpos>
            map_rs_to_chrpos( file_reading:: GenotypeFileHandle raw_ref_file);
static
std:: unordered_multimap<chrpos, string>
            make_map_chrpos_to_SNPname( file_reading:: GenotypeFileHandle raw_ref_file);
static
vector<vector<int>>
lookup_many_ref_rows( vector<SNPiterator<GenotypeFileHandle>>   const &  snps
                , file_reading:: CacheOfRefPanelData              &  cache
                );
static
mvn:: SquareMatrix
make_C_tag_tag_matrix(
                    vector<vector<int>>              const & genotypes_for_the_tags
                    , double lambda
                    );
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int>>              const & genotypes_for_the_tags
        , vector<vector<int>>              const & genotypes_for_the_unks
        , vector<SNPiterator<GenotypeFileHandle>> const & tag_its
        , vector<SNPiterator<GenotypeFileHandle>> const & unk_its
        , double                                          lambda
        );

enum class which_direction_t { DIRECTION_SHOULD_BE_REVERSED
                             , NO_ALLELE_MATCH
                             , DIRECTION_AS_IS };
string to_string(which_direction_t dir) {
    switch(dir) {
        break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED: return "DIRECTION_SHOULD_BE_REVERSED";
        break; case which_direction_t:: NO_ALLELE_MATCH:              return "NO_ALLELE_MATCH";
        break; case which_direction_t:: DIRECTION_AS_IS:              return "DIRECTION_AS_IS";
    }
    assert(1==2);
    return ""; // should never get here
}
static
    which_direction_t decide_on_a_direction
    ( SNPiterator<GenotypeFileHandle> const & //r
    , SNPiterator<GwasFileHandle>     const & //g
    );
} // namespace ssimp

static
void    set_appropriate_locale(ostream & stream) {
    // By default, apply the user's locale, for example the thousands separator ...
    stream.imbue(std::locale(""));

    /*  However, within our test scripts, we force apostrophe ['] as the thousands
     *  separator so that we can compare results consistently across different developers
     *
     *  More info: http://en.cppreference.com/w/cpp/locale/numpunct/thousands_sep
     */

    if(std:: getenv("FORCE_THOUSANDS_SEPARATOR") != nullptr) {
        auto force_thousands_separator = std::string(std:: getenv("FORCE_THOUSANDS_SEPARATOR"));
        if(force_thousands_separator.size() ==1) {
            struct custom_thousands_separator : std::numpunct<char> {
                char    m_custom_separator;

                custom_thousands_separator(char c) : m_custom_separator(c) {}

                char do_thousands_sep()   const { return m_custom_separator; }
                std::string do_grouping() const { return "\3"; } // groups of 3 digits
            };
            stream.imbue(   std::locale (   stream.getloc()
                                        ,   new custom_thousands_separator(force_thousands_separator.at(0))
                                        // Nobody likes 'new' nowadays, but apparently this won't
                                        // leak. From http://en.cppreference.com/w/cpp/locale/locale/locale:
                                        //
                                        //  "Overload 7 is typically called with its second argument, f, obtained directly from a
                                        //   new-expression: the locale is responsible for calling the matching delete from its
                                        //   own destructor."
                                        )
                        );
        }
    }
}

int main(int argc, char **argv) {

    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);

    set_appropriate_locale(cout);

    cout << std::setprecision(20);

    if( options:: opt_raw_ref.empty() ||  options:: opt_gwas_filename.empty()) {
        DIE("Should pass args.\n    Usage:   " + string(argv[0]) + " --ref REFERENCEVCF --gwas GWAS --lambda 0.0 --window.width 1000000 --flanking.width 250000");
    }

    if( !options:: opt_impute_snps.empty() ) {
        gz:: igzstream f(options:: opt_impute_snps.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << options:: opt_impute_snps << ']');
        string one_SNPname_or_chrpos;
        while(getline( f, one_SNPname_or_chrpos)) {
            options:: opt_impute_snps_as_a_uset.insert(one_SNPname_or_chrpos);
        }
        options:: opt_impute_snps_as_a_uset.empty() && DIE("Nothing in --impute.snps file? [" << options:: opt_impute_snps << "]");
    }

    if(!options:: opt_raw_ref.empty() && !options:: opt_gwas_filename.empty()) {
        PP( options:: opt_raw_ref
          , options:: opt_gwas_filename
          , options:: opt_window_width
          );

        // Load the two files
        auto raw_ref_file = file_reading:: read_in_a_raw_ref_file(options:: opt_raw_ref);
        auto gwas         = file_reading:: read_in_a_gwas_file(options:: opt_gwas_filename);

        // Compare the chrpos in both files.
        // If GWAS has chrpos, it should be the same as in the RefPanel.
        // If GWAS doesn't have chrpos, copy it from the ref panel
        auto m            = ssimp:: map_rs_to_chrpos( raw_ref_file );
        update_positions_by_comparing_to_another_set( gwas, m );

        // The RefPanel has already been automatically sorted by chrpos,
        // but now we must do it for the GWAS, as the positions might have
        // changed in the GWAS.
        gwas->sort_my_entries();

        // Count how many GWAS SNPs have no position, and therefore
        // will be ignored.
        auto number_of_GWASsnps_with_unknown_position = std:: count_if(
                begin_from_file(gwas)
               ,  end_from_file(gwas)
               , [](auto v) { return v == chrpos{-1,-1}; }
               );

        // Print the various SNP counts
        cout << '\n';
        PP(raw_ref_file->number_of_snps());
        PP(        gwas->number_of_snps());
        PP(number_of_GWASsnps_with_unknown_position);

        cout << '\n';
        // Go through regions, printing how many
        // SNPs there are in each region
        ssimp:: impute_all_the_regions(raw_ref_file, gwas);
    }
}

namespace ssimp{

using file_reading:: GenotypeFileHandle;
using file_reading:: GwasFileHandle;

static
bool decide_whether_to_skip_this_tag( SNPiterator<GenotypeFileHandle>       const & it ) {
    if( options:: opt_impute_range.empty()
     && options:: opt_impute_snps.empty()
     )
        return false;

    // You can't apply both. Only one or the other:
    assert( options:: opt_impute_range.empty() != options:: opt_impute_snps.empty());

    if(!options:: opt_impute_snps.empty()) {
        // --impute.snps was specified, so the file has been loaded into ...
        assert(!options:: opt_impute_snps_as_a_uset.empty());
        return  (   options:: opt_impute_snps_as_a_uset.count( it.get_SNPname() )
                +   options:: opt_impute_snps_as_a_uset.count( AMD_FORMATTED_STRING("chr{0}:{1}", it.get_chrpos().chr, it.get_chrpos().pos ) )
                ) == 0;
    }

    if(!options:: opt_impute_range.empty()) {
        // impute.range will either be an entire chromosome,
        // or two locations seperated by a '-'

        auto separated_by_hyphen = utils:: tokenize(options:: opt_impute_range, '-');

        // remove any leading 'chr' or 'Chr'
        for(auto &s : separated_by_hyphen) {
            if(     s.substr(0,3) == "chr"
                 || s.substr(0,3) == "Chr")
                s = s.substr(3);
        };

        auto lambda_chrpos_text_to_object = [](string const &as_text, bool to_end_of_chromosome) {
            auto separated_by_colon = utils:: tokenize(as_text,':');
            chrpos position;
            switch(separated_by_colon.size()) {
                break; case 1: // no colon, just a chromosome
                        position.chr = utils:: lexical_cast<int>(separated_by_colon.at(0));
                        position.pos = to_end_of_chromosome
                                        ?  std:: numeric_limits<int>::max()
                                        :  std:: numeric_limits<int>::lowest() ;
                break; case 2:
                        position.chr = utils:: lexical_cast<int>(separated_by_colon.at(0));
                        position.pos = utils:: lexical_cast<int>(separated_by_colon.at(1));
                break; default:
                        DIE("too many colons in [" << as_text << "]");
            }
            return position;
        };

        switch(separated_by_hyphen.size()) {
            break; case 1:
            {
                chrpos  lower_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), false);
                chrpos  upper_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), true);
                return  !(  it.get_chrpos() >= lower_allowed
                         && it.get_chrpos() <= upper_allowed    );
            }
            break; case 2:
            {
                chrpos  lower_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), false);
                chrpos  upper_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(1).c_str(), true);
                return  !(  it.get_chrpos() >= lower_allowed
                         && it.get_chrpos() <= upper_allowed    );
            }
           break; default:
                        DIE("too many hyphens in [" << options:: opt_impute_range << "]");
        }
    }

    assert(0); // shouldn't get here
    return false;
}

static
void impute_all_the_regions( file_reading:: GenotypeFileHandle         ref_panel
                             , file_reading:: GwasFileHandle             gwas
                             ) {
    unique_ptr<ofstream>   out_stream_for_imputations;
    ostream              * out_stream_ptr = nullptr;
    if(!options:: opt_out.empty()) {
        out_stream_for_imputations = make_unique<ofstream>( options:: opt_out. c_str() );
        (*out_stream_for_imputations) || DIE("Couldn't open the --out file [" + options:: opt_out + "]");
        out_stream_ptr = &*out_stream_for_imputations;
    }
    else {
        // If --out isn't specified, just print to stdout
        out_stream_ptr = & cout;
    }
    assert(out_stream_ptr);
    {   // print the header line
                    (*out_stream_ptr)
                        << "chr:pos"
                        << '\t' << "z_imp"
                        << '\t' << "SNPname"
                        << '\t' << gwas->get_column_name_allele_ref() // copy column name from the GWAS input
                        << '\t' << gwas->get_column_name_allele_alt() // copy column name from the GWAS input
                        << '\t' << "impqual"
                        << endl;
    }

    auto const b_ref  = begin_from_file(ref_panel);
    auto const e_ref  =   end_from_file(ref_panel);
    auto const b_gwas = begin_from_file(gwas);
    auto const e_gwas =   end_from_file(gwas);

    PP(options:: opt_window_width
      ,options:: opt_flanking_width
            );

    auto const map_chrpos_to_SNPname           = ssimp:: make_map_chrpos_to_SNPname( ref_panel );

    tbi:: read_vcf_with_tbi ref_vcf {
        "ref/tbi/TWINSUK.every100thRS.chrm123.100people.vcf.gz"
        //options:: opt_raw_ref
    };

    for(int chrm =  1; chrm <= 22; ++chrm) {
        file_reading:: CacheOfRefPanelData cache(ref_panel);;

        // First, find the begin and end of this chromosome
        auto c_begin = std:: lower_bound(b_ref, e_ref, chrpos{chrm, std::numeric_limits<int>::lowest() });
        auto c_end   = std:: lower_bound(b_ref, e_ref, chrpos{chrm, std::numeric_limits<int>::max()  });
        assert(c_end >= c_begin);
        if(c_begin != c_end)
            assert(c_begin.get_chrpos().pos >= 0); // first position is at least zero

        for(int w = 0; ; ++w ) {
            int current_window_start = w     * options:: opt_window_width;
            int current_window_end   = (w+1) * options:: opt_window_width;

            /*
             * For a given window, there are three "ranges" to consider:
             *  -   The range of GWAS SNPs in the "broad" window (i.e. including the flanking region)
             *  -   The range of reference SNPs in the "broad" window.
             *  -   The range of rererence SNPs in the "narrow" window (i.e. without the flanking region)
             *
             *  The intersection of the first two of these three ranges is the set of SNPs that
             *  are useful as tags. The third range is the set of targets.
             *
             *  (For now, we're including every tag as a target too.
             *
             *  In the next few lines, we find the endpoints (begin and end) of each of these
             *  three ranges.
             */

            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            auto w_ref_wide_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_ref_wide_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end   + options:: opt_flanking_width });

            auto w_ref_narrow_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start});
            auto w_ref_narrow_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end  });

            struct RefRecord {
                int     pos;
                string  ID;
                string  ref;
                string  alt;
                vector<int> z12;

                bool operator< (file_reading:: chrpos crps) {
                                                                    return pos < crps.pos;
                                                                }
            };

            vector<RefRecord>   all_nearby_ref_data;
            {   // Read in all the reference panel data in this broad window, but bi-allele SNPs.
                ref_vcf.set_region  (   chrpos{chrm,current_window_start - options:: opt_flanking_width}
                                    ,   chrpos{chrm,current_window_end   + options:: opt_flanking_width}
                                    );
                VcfRecord record;
                int N = -1; // just to verify it is constant
                while(ref_vcf.reader.readRecord(record)) {
                    RefRecord rr;
                    rr.pos      =   record.get1BasedPosition();
                    rr.ID       =   record.getIDStr();
                    rr.ref      =   record.getRefStr();
                    rr.alt      =   record.getAltStr();
                    assert(record.getNumAlts() == 1); // The 'DISCARD' rule should already have skipped those with more alts
                    assert(record.hasAllGenotypeAlleles());
                    if(N==-1) {
                        N = record.getNumSamples();
                    }
                    assert(N == record.getNumSamples());
                    for(int i=0;i<N;++i) {
                        assert(2==record.getNumGTs(i));
                        int l = record.getGT(i, 0);
                        int r = record.getGT(i, 1);
                        assert((l | 1) == 1);
                        assert((r | 1) == 1);
                        rr.z12.push_back(l+r);
                    }
                    assert(N == ssize(rr.z12));

                    all_nearby_ref_data.push_back(rr);

                    if(0)
                    PP( record.getChromStr()
                      , record.get1BasedPosition()
                      , record.getIDStr()
                      , record.getNumSamples()         // TODO: Check this is constant?
                      , record.hasAllGenotypeAlleles() // TODO: Should always be true?
                      );

                }

                for(int o=0; o+1 < ssize(all_nearby_ref_data); ++o) { // check position is monotonic
                    assert(all_nearby_ref_data.at(o).pos <= all_nearby_ref_data.at(o+1).pos);
                    assert(all_nearby_ref_data.at(o).pos <  all_nearby_ref_data.at(o+1).pos); // TODO: remove this, after making sure it breaks somewhere!
                }
            }

            // Check if the current target window is past the end of the chromosome, in which case
            // we 'break' out in order to allow it to finish this chromosome and move onto
            // the next chromosome
            if(w_ref_narrow_begin == c_end)
                break; // Finished with this chromosome

            // If the current target window has no SNPs, then 'continue' to the next window
            if(w_ref_narrow_begin == w_ref_narrow_end)
                continue; // Nothing to impute, skip to the next window


            // Next, look up each 'tag candidate' in turn and - if it's a suitable SNP - store
            // it's z-statistic and an iterator to the reference panel data for the same SNP.
            vector<double>                          tag_zs;
            vector<SNPiterator<GenotypeFileHandle>> tag_its;
            for(auto tag_candidate = range:: range_from_begin_end( w_gwas_begin, w_gwas_end )
                    ; ! tag_candidate.empty()
                    ;   tag_candidate.advance()
                    ) {
                auto crps = tag_candidate.current_it().get_chrpos();
                // Find the ref panel entries at the exact same chr:pos
                auto ref_candidates = range:: range_from_begin_end(
                         std:: lower_bound( w_ref_wide_begin , w_ref_wide_end, crps )
                        ,std:: upper_bound( w_ref_wide_begin , w_ref_wide_end, crps )
                        );

                // There may be multiple reference panel SNPs at the same chr:pos, for
                // example an 'indel' as well as a conventional biallelic SNP.
                // These next few lines find the collection of reference panel SNPs
                // that have appropriate alleles.
                vector<double> one_tag_zs;
                vector<SNPiterator<GenotypeFileHandle>> one_tag_its;
                for(; ! ref_candidates.empty(); ref_candidates.advance()) {
                    assert( ref_candidates.current_it().get_chrpos() == crps );

                    // If there is more than one alt allele, we need to skip this,
                    // *and* if there is no variation on the 012
                    if(cache.lookup_one_ref_get_calls(ref_candidates.current_it()) .size() == 0) {
                        continue;
                    }

                    auto dir = decide_on_a_direction( ref_candidates.current_it()
                                         , tag_candidate .current_it() );
                    switch(dir) {
                        break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED:
                            one_tag_zs.push_back( -tag_candidate.current_it().get_z() );
                            one_tag_its.push_back( ref_candidates.current_it()        );
                        break; case which_direction_t:: DIRECTION_AS_IS             :
                            one_tag_zs.push_back(  tag_candidate.current_it().get_z() );
                            one_tag_its.push_back( ref_candidates.current_it()        );
                        break; case which_direction_t:: NO_ALLELE_MATCH             : ;
                    }
                }
                // Finally, we make the decision on whether to include this tag or not.
                // A tag is useful if there is exactly one reference SNP that corresponds
                // to it (where 'correspond' means that the positions and alleles match).
                assert(one_tag_zs.size() == one_tag_its.size());
                if(1==one_tag_zs.size()) {
                    for(auto z : one_tag_zs)
                        tag_zs.push_back(z);
                    for(auto it : one_tag_its)
                        tag_its.push_back(it);
                }
            }
            assert(tag_zs.size() == tag_its.size());

            int const number_of_tags = tag_zs.size();
            if(number_of_tags == 0)
                continue;

            // Now, find suitable targets - i.e. anything in the reference panel in the narrow window
            // But some SNPs will have to be dropped, depending on --impute.range and --impute.snps
            vector<SNPiterator<GenotypeFileHandle>> unk_its;
            vector<RefRecord*>                      unk2_its;
            for(auto it = w_ref_narrow_begin; it<w_ref_narrow_end; ++it) {
                // actually, we should think about ignoring SNPs in certain situations

                bool skip_this = decide_whether_to_skip_this_tag(it);
                if(skip_this)
                    continue;

                auto const & z12_for_this_SNP = cache.lookup_one_ref_get_calls(it);
                if (z12_for_this_SNP.empty()) {
                    // Empty vector means the SNP is not binary in the reference panel.
                    // For now, we're ignoring such SNPs, but maybe we could include them
                    // by simply using 'NA' in place of the awkward individuals.
                    continue;
                }
                auto z12_minmax = minmax_element(z12_for_this_SNP.begin(), z12_for_this_SNP.end());
                if  (*z12_minmax.first == *z12_minmax.second){
                    continue;   // no variation in the allele counts for this SNP within
                                //the ref panel, therefore useless for imputation
                }
                unk_its         .push_back( it              );
            }
            {
                auto w_ref_narrow_begin = std:: lower_bound (   all_nearby_ref_data.begin()
                                                            ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_start}
                                                            //,   [](RefRecord &l, file_reading:: chrpos crps) { return l.pos < crps.pos; }
                                                            );
                auto w_ref_narrow_end = std:: lower_bound (   all_nearby_ref_data.begin()
                                                            ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_end}
                                                            //,   [](RefRecord &l, file_reading:: chrpos crps) { return l.pos < crps.pos; }
                                                            );
                for(auto it = w_ref_narrow_begin; it<w_ref_narrow_end; ++it) {
                    // actually, we should think about ignoring SNPs in certain situations

                    //bool skip_this = decide_whether_to_skip_this_tag(it); if(skip_this) continue; // TODO: implement this

                    auto const & z12_for_this_SNP = it->z12;
                    auto z12_minmax = minmax_element(z12_for_this_SNP.begin(), z12_for_this_SNP.end());
                    if  (*z12_minmax.first == *z12_minmax.second){
                        continue;   // no variation in the allele counts for this SNP within
                                    //the ref panel, therefore useless for imputation
                    }
                    unk2_its    .push_back( &*it       );
                }
            }
#if 0
            for(auto && l : unk_its) {
                PP(l.get_chrpos().pos
                  ,l.get_SNPname()
                  ,l.get_allele_ref()
                  ,l.get_allele_alt()   );
            }
            for(auto && r : unk2_its) {
                PP(r->pos, r->ID, r->ref, r->alt);
            }
            PP(unk_its.size() , unk2_its.size());
            for(auto && x : range:: zip( from::vector(unk_its), from::vector(unk2_its))) {
                auto && l = x | stdget0;
                auto && r = x | stdget1;
                PP    (l.get_chrpos().pos   ,       r->pos);
                assert(l.get_chrpos().pos   ==      r->pos);
                assert(l.get_SNPname()      ==      r->ID);
                assert(l.get_allele_ref()   ==      r->ref);
                assert(l.get_allele_alt()   ==      r->alt);
            }
            PP(unk_its.size() , unk2_its.size());
#endif
            assert(unk_its.size() == unk2_its.size());


            // Next few lines are kind of boring, just printing some statistics.

            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;
            int  number_of_all_targets                          = unk_its.size();

            if(number_of_all_targets == 0)
                continue;

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << number_of_tags                                 << " # SNPs in both (i.e. useful as tags)\n";
            cout << setw(8) << number_of_all_targets                          << " # target SNPs (anything in narrow window, will include some tags)\n";

            // Next, actually look up all the relevant SNPs within the reference panel
            vector<vector<int>> genotypes_for_the_tags = lookup_many_ref_rows(tag_its, cache);
            vector<vector<int>> genotypes_for_the_unks = lookup_many_ref_rows(unk_its, cache);
            vector<vector<int>> genotypes_for_the_unks_ = from:: vector(unk2_its) |view::map| [](RefRecord *rrp) { return rrp->z12; } |action::collect;
            assert(genotypes_for_the_unks == genotypes_for_the_unks_);

            assert(number_of_tags        == utils:: ssize(genotypes_for_the_tags));
            assert(number_of_all_targets == utils:: ssize(genotypes_for_the_unks));

            static_assert( std:: is_same< vector<vector<int>> , decltype(genotypes_for_the_tags) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation
            static_assert( std:: is_same< vector<vector<int>> , decltype(genotypes_for_the_unks) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation

            assert(!genotypes_for_the_tags.empty());
            assert(!genotypes_for_the_unks.empty());

            int const N_ref = genotypes_for_the_tags.at(0).size(); // the number of individuals
            assert(N_ref > 0);


            /*
             * Next few lines do a lot. The compute correlation, applying lambda regularization, and do imputation:
             */
            mvn:: SquareMatrix  C           = make_C_tag_tag_matrix(genotypes_for_the_tags, options:: opt_lambda);
            mvn:: VecCol        C_inv_zs    = solve_a_matrix (C, mvn:: make_VecCol(tag_zs));
            mvn:: Matrix        c           = make_c_unkn_tags_matrix( genotypes_for_the_tags
                                                         , genotypes_for_the_unks
                                                         , tag_its
                                                         , unk_its
                                                         , options:: opt_lambda
                                                         );

            // Next two lines, are the imputation, and its quality
            auto c_Cinv_zs = mvn:: multiply_matrix_by_colvec_giving_colvec(c, C_inv_zs);
            auto   Cinv_c  =     mvn:: multiply_NoTrans_Trans(invert_a_matrix (C) , c);

            assert( (int)Cinv_c.size1() == number_of_tags );
            assert( (int)Cinv_c.size2() == number_of_all_targets );


            // Finally, print out the imputations
            assert(number_of_all_targets == ssize(unk_its));
            mvn:: Matrix to_store_one_imputation_quality(1,1); // a one-by-one-matrix
            for(int i=0; i<number_of_all_targets; ++i) {
                    auto && target = unk_its.at(i);
                    auto crps = target.get_chrpos();
                    auto SNPname = target.get_SNPname();

                    auto imp_qual = [&](){ // multiply one vector by another, to directly get
                                            // the diagonal of the imputation quality matrix.
                                            // This should be quicker than fully computing
                                            // c' * inv(C) * c
                        auto lhs = gsl_matrix_const_submatrix(     c.get(), i, 0, 1, number_of_tags);
                        auto rhs = gsl_matrix_const_submatrix(Cinv_c.get(), 0, i, number_of_tags, 1);
                        const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &lhs.matrix, &rhs.matrix, 0, to_store_one_imputation_quality.get());
                        assert(res_0 == 0);
                        return to_store_one_imputation_quality(0,0);
                    }();

                    (*out_stream_ptr)
                        << crps
                        << '\t' << c_Cinv_zs(i)
                        << '\t' << SNPname
                        << '\t' << target.get_allele_ref()
                        << '\t' << target.get_allele_alt()
                        << '\t' << imp_qual
                        << endl;
            }
        }
    }
}

static
std:: unordered_map<string, chrpos>
            map_rs_to_chrpos( file_reading:: GenotypeFileHandle raw_ref_file ) {
    auto       b = begin_from_file(raw_ref_file);
    auto const e =   end_from_file(raw_ref_file);
    std:: unordered_map<string, chrpos> m;
    for(;b<e; ++b) {
        if(b.get_SNPname() == ".") // skip these
            continue;
        auto rel = m.insert( std:: make_pair(b.get_SNPname(), b.get_chrpos()) );
        rel.second || DIE("same SNPname twice in the ref panel [" << b.get_SNPname() << "]");
    }
    return m;
}
static
std:: unordered_multimap<chrpos, string>
            make_map_chrpos_to_SNPname( file_reading:: GenotypeFileHandle raw_ref_file ) {
    auto       b = begin_from_file(raw_ref_file);
    auto const e =   end_from_file(raw_ref_file);
    std:: unordered_multimap<chrpos, string> m;
    for(;b<e; ++b) {
        m.insert( std:: make_pair(b.get_chrpos(), b.get_SNPname()) );
    }
    return m;
}
static
vector<vector<int>>
lookup_many_ref_rows( vector<SNPiterator<GenotypeFileHandle>>   const &  snps
                , file_reading:: CacheOfRefPanelData              &  cache
                ) {
    vector<vector<int>> many_calls;
    for(auto crps : snps) {
        many_calls.push_back( cache.lookup_one_ref_get_calls(crps));
    }
    assert(many_calls.size() == snps.size());
    return many_calls;
}

static
mvn:: SquareMatrix
make_C_tag_tag_matrix( vector<vector<int>>              const & genotypes_for_the_tags
                     , double                                   lambda
                     ) {
    int const number_of_tags = genotypes_for_the_tags.size();
    int const N_ref = genotypes_for_the_tags.at(0).size();
    assert(N_ref > 0);

    mvn:: SquareMatrix C (number_of_tags);
    for(int k=0; k<number_of_tags; ++k) {
        for(int l=0; l<number_of_tags; ++l) {
            assert(N_ref        == utils:: ssize(genotypes_for_the_tags.at(k)));
            assert(N_ref        == utils:: ssize(genotypes_for_the_tags.at(l)));
            double c_kl = gsl_stats_int_correlation( &genotypes_for_the_tags.at(k).front(), 1
                                                   , &genotypes_for_the_tags.at(l).front(), 1
                                                   , N_ref );
            if(c_kl > (1.0-1e-10)) { // sometimes it sneaks above one, don't really know how
                assert(c_kl-1.0 < 1e-5);
                assert(c_kl-1.0 >-1e-5);
                c_kl = 1.0;
            }
            if(c_kl < (-1.0+1e-10)) { // in case it goes below -1.0 also
                assert(c_kl+1.0 < 1e-5);
                assert(c_kl+1.0 >-1e-5);
                c_kl =-1.0;
            }
            assert(c_kl >= -1.0);
            assert(c_kl <=  1.0);
            if(k==l) {
                assert(c_kl == 1.0);
                C.set(k,l,c_kl);
            }
            else {
                C.set(k,l,c_kl * (1.0-lambda));
            }
        }
    }
    return C;
}
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int>>              const & genotypes_for_the_tags
        , vector<vector<int>>              const & genotypes_for_the_unks
        , vector<SNPiterator<GenotypeFileHandle>> const & tag_its
        , vector<SNPiterator<GenotypeFileHandle>> const & unk_its
        , double                                          lambda
        ) {
    int const number_of_tags = genotypes_for_the_tags.size();
    int const number_of_all_targets = genotypes_for_the_unks.size();
    assert(number_of_tags        == ssize(tag_its));
    assert(number_of_all_targets == ssize(unk_its));
    int const N_ref = genotypes_for_the_tags.at(0).size();
    assert(N_ref > 0);

    mvn:: Matrix      c(number_of_all_targets, number_of_tags);

    for(auto tags : zip_val ( range:: ints(number_of_tags)
                            , range:: range_from_begin_end(tag_its)
                            , range:: range_from_begin_end(genotypes_for_the_tags) | view:: ref_wraps
                )) {
        int     k           = std::get<0>(tags);
        auto    tag_its_k   = std::get<1>(tags);
        auto &  calls_at_k  = std::get<2>(tags) .get();

        assert(&calls_at_k == &genotypes_for_the_tags.at(k)); // double check that we're getting it by-reference, i.e via the .get() on the previous line

        assert(tag_its_k == tag_its.at(k));

        zip_val ( range:: ints(number_of_all_targets)
                , range:: range_from_begin_end(unk_its)
                , range:: range_from_begin_end(genotypes_for_the_unks) | view:: ref_wraps
                )
        |action:: unzip_foreach|
        [&] (   int
                    u
            ,   SNPiterator<GenotypeFileHandle>
                    unk_its_u
            ,   std:: reference_wrapper< vector<int> const >
                    calls_at_u_rw
            ) -> void {
                vector<int> const &  calls_at_u  = calls_at_u_rw .get();
                assert(&calls_at_u == &genotypes_for_the_unks.at(u));

                assert(N_ref        == utils:: ssize(calls_at_k));
                assert(N_ref        == utils:: ssize(calls_at_u));
                double c_ku = gsl_stats_int_correlation( &calls_at_k.at(0), 1
                                                       , &calls_at_u.at(0), 1
                                                       , N_ref );
                if(c_ku > (1.0-1e-10)) {
                    assert(c_ku-1.0 < 1e-5);
                    assert(c_ku-1.0 >-1e-5);
                    c_ku = 1.0;
                }
                if(c_ku < (-1.0+1e-10)) {
                    assert(c_ku+1.0 < 1e-5);
                    assert(c_ku+1.0 >-1e-5);
                    c_ku = -1.0;
                }

                if(tag_its_k.m_line_number == unk_its_u.m_line_number) {
                    assert(c_ku ==  1.0);
                    c.set(u,k,c_ku);
                }
                else {
                    assert(c_ku >= -1.0);
                    assert(c_ku <=  1.0);
                    c.set(u,k,c_ku * (1.0-lambda));
                }
         };
    }
    return c;
};
static
    which_direction_t decide_on_a_direction
    ( SNPiterator<GenotypeFileHandle> const & r
    , SNPiterator<GwasFileHandle>     const & g
    ) {
        auto const rp_ref = r.get_allele_ref();
        auto const rp_alt = r.get_allele_alt();
        auto const gw_ref = g.get_allele_ref();
        auto const gw_alt = g.get_allele_alt();

        //PP(r.get_chrpos() ,rp_ref ,rp_alt ,gw_ref ,gw_alt);
        assert(rp_ref != rp_alt);
        assert(gw_ref != gw_alt);

        if( rp_ref == gw_ref
         && rp_alt == gw_alt)
            return which_direction_t:: DIRECTION_AS_IS;

        if( rp_ref == gw_alt
         && rp_alt == gw_ref)
            return which_direction_t:: DIRECTION_SHOULD_BE_REVERSED;

        return which_direction_t:: NO_ALLELE_MATCH;
    }

} // namespace ssimp
