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

using utils:: ssize;
using utils:: print_type;
using utils:: stdget0;
using utils:: stdget1;

using tbi:: RefRecord;


namespace ssimp {
// Some forward declarations

static
void impute_all_the_regions(   string                                   filename_of_vcf
                             , file_reading:: GwasFileHandle_NONCONST   gwas
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
        , vector<RefRecord const *              > const & tag_its
        , vector<RefRecord const *              > const & unk_its
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
    ( RefRecord                       const & //r
    , SNPiterator     const & //g
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
        auto gwas         = file_reading:: read_in_a_gwas_file(options:: opt_gwas_filename);

        // Go through regions, printing how many
        // SNPs there are in each region
        ssimp:: impute_all_the_regions(options:: opt_raw_ref, gwas);
    }
}

namespace ssimp{

struct skipper_target_I {
    virtual bool    skip_me(int chrm, RefRecord const * rrp)    =0;
};

static
chrpos lambda_chrpos_text_to_object     (string const &as_text, bool to_end_of_chromosome) {
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
}

unique_ptr<skipper_target_I> make_skipper_for_targets
        (   std:: string                                    range_as_string
        ,   std::unordered_set<std::string>     const &     uset_of_strings
        ) {
    if( range_as_string.empty() && uset_of_strings.empty())
    {
        struct local : skipper_target_I {
            bool    skip_me(int     , RefRecord const *    )   override { return false; }
        };
        return make_unique<local>();
    }
    if( range_as_string.empty() && !uset_of_strings.empty())
    {   // --impute.snps was specified
        // We check if either the ID, or the chr:pos is in the --impute.snps file
        struct local : skipper_target_I {
            decltype(uset_of_strings) const & m_uset_of_strings;
            local(decltype(m_uset_of_strings) const & r) : m_uset_of_strings(r) {}

            bool    skip_me(int chrm, RefRecord const * rrp)   override {
                assert(!m_uset_of_strings.empty());
                return  (   m_uset_of_strings.count( rrp->ID )
                        +   m_uset_of_strings.count( AMD_FORMATTED_STRING("chr{0}:{1}", chrm, rrp->pos ) )
                        ) == 0;
            }
        };
        return make_unique<local>(uset_of_strings);
    }
    if( !range_as_string.empty() && uset_of_strings.empty())
    {   // --impute.range was specified
        // --impute.range will either be an entire chromosome,
        // or two locations in one chromosome separated by a '-'

        auto separated_by_hyphen = utils:: tokenize(range_as_string, '-');

        // remove any leading 'chr' or 'Chr'
        for(auto &s : separated_by_hyphen) {
            if(     s.substr(0,3) == "chr"
                 || s.substr(0,3) == "Chr")
                s = s.substr(3);
        };

        struct local : skipper_target_I {
            chrpos  lower_allowed;
            chrpos  upper_allowed;
            bool    skip_me(int chrm, RefRecord const * rrp)   override {
                bool b = !(  chrpos{chrm, rrp->pos} >= lower_allowed && chrpos{chrm, rrp->pos} <= upper_allowed    );
                return b;
            }
        };
        auto up = make_unique<local>();

        switch(separated_by_hyphen.size()) {
            break; case 1:
            {
                up->lower_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), false);
                up->upper_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), true);
            }
            break; case 2:
            {
                up->lower_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(0).c_str(), false);
                up->upper_allowed = lambda_chrpos_text_to_object(separated_by_hyphen.at(1).c_str(), true);
            }
           break; default:
                        DIE("too many hyphens in [" << range_as_string << "]");
        }
        return std::move(up);
    }

    assert  (   range_as_string.empty()
             && uset_of_strings.empty());
    DIE("--impute.range and --impute.snps can't be used together. (I thought I had already checked for this)");
    return nullptr;
    // TODO: clean up the end of this function
    assert(0);
        struct local : skipper_target_I {
            bool    skip_me(int     , RefRecord const *    )   override { return false; }
        };
        return make_unique<local>();
}

static
void impute_all_the_regions(   string                                   filename_of_vcf
                             , file_reading:: GwasFileHandle_NONCONST   gwas
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


    PP(options:: opt_window_width
      ,options:: opt_flanking_width
            );

    auto skipper_target = make_skipper_for_targets(options:: opt_impute_range, options:: opt_impute_snps_as_a_uset);

    tbi:: read_vcf_with_tbi ref_vcf { filename_of_vcf };

    for(int chrm =  1; chrm <= 22; ++chrm) {

        for(int w = 0; ; ++w ) {
            int current_window_start = w     * options:: opt_window_width;
            int current_window_end   = (w+1) * options:: opt_window_width;

            /*
             * For a given window, there are three "ranges" to consider:
             *  -   The range of GWAS SNPs in the "broad" window (i.e. including the flanking region)
             *  -   The range of reference SNPs in the "broad" window.
             *  -   The range of reference SNPs in the "narrow" window (i.e. without the flanking region)
             *
             *  The intersection of the first two of these three ranges is the set of SNPs that
             *  are useful as tags. The third range is the set of targets.
             *
             *  (For now, we're including every tag as a target too.
             *
             *  In the next few lines, we find the endpoints (begin and end) of each of these
             *  three ranges.
             */

            // First, load up all the reference data in the broad window.
            // If this is empty, and it's the last window, then we can finish with this chromosome.

            vector<RefRecord>   all_nearby_ref_data;
            bool not_the_last_window = false;
            {   // Read in all the reference panel data in this broad window, only bi-allelic SNPs.
                ref_vcf.set_region  (   chrpos{chrm,current_window_start - options:: opt_flanking_width}
                                    ,   chrpos{chrm,std::numeric_limits<int>::max()                    } // in theory, to the end of the chromosome. See a few lines below
                                    );
                RefRecord rr;
                while(ref_vcf.read_record_into_a_RefRecord(rr)) {

                    /*
                     * if the position is beyond the current wide window, it means
                     * this is *not* the last window
                     */
                    if(rr.pos   >= current_window_end   + options:: opt_flanking_width) {
                        not_the_last_window = true;
                        break;
                    }

                    all_nearby_ref_data.push_back(rr);

                }

                for(int o=0; o+1 < ssize(all_nearby_ref_data); ++o) { // check position is monotonic
                    assert(all_nearby_ref_data.at(o).pos <= all_nearby_ref_data.at(o+1).pos);
                }
            }

            if(all_nearby_ref_data.empty() && !not_the_last_window) {
                // No more reference panel data, therefore no tags, therefore end of chromosome
                break;
            }

            if(all_nearby_ref_data.empty()) {
                // no tags
                continue; // to the next window in this chromosome
            }

            {   // Update the gwas data with position information from the ref panel
                // I should put this chunk of code into another function
                std:: unordered_map<string, int> map_SNPname_to_pos; // just for this broad window
                for(auto & rr : all_nearby_ref_data) {
                    switch(map_SNPname_to_pos.count(rr.ID)) {
                        break; case 0:
                            map_SNPname_to_pos[rr.ID] = rr.pos;
                        break; case 1:
                        {
                            if(map_SNPname_to_pos[rr.ID] != rr.pos) {
                                PP(rr.ID); // TODO: remove this. I think it will be ID='.'
                                map_SNPname_to_pos[rr.ID] = -1;
                            }
                        }
                    }
                }
                for(auto & msp : map_SNPname_to_pos) {
                    assert(msp.second != -1); // TODO: remove these instead
                }

                // Now that we have the position data from ref panel, see if we can copy it into the GWAS data
                auto gwas_it = begin_from_file(gwas);
                auto const e_gwas =   end_from_file(gwas);
                for(; gwas_it != e_gwas; ++gwas_it) {
                    auto gwas_SNPname = gwas_it.get_SNPname();
                    if(                     map_SNPname_to_pos.count( gwas_SNPname )) {
                        auto pos_in_ref =   map_SNPname_to_pos      [ gwas_SNPname ];
                        if  (   gwas_it.get_chrpos().chr == chrm
                             && gwas_it.get_chrpos().pos == pos_in_ref) {
                            // this fine, nothing to change, the two data sources agree
                        }
                        else if (   gwas_it.get_chrpos().chr == -1
                                 && gwas_it.get_chrpos().pos == -1) {
                            // just copy it in
                                gwas_it.set_chrpos(chrpos{chrm,pos_in_ref});
                        }
                        else
                            DIE("disagreement on position");
                    }
                }
            }
            // Now, thanks to the block just completed, we have all the positions in the GWAS
            // that are relevant for this window.
            gwas->sort_my_entries(); // Sort them by position

            // Of the next four lines, only the two are interesting (to find the range of GWAS snps for this window)
            auto const b_gwas = begin_from_file(gwas);
            auto const e_gwas =   end_from_file(gwas);
            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            // Next, find tags
            vector<double>                          tag_zs_;
            vector<RefRecord const *              > tag_its_;
            {   // Look for tags in the broad window.
                // This means moving monotonically through 'all_nearby_ref_data' and
                // the gwas data in parallel, and identifying suitable 'pairs'
                for(auto tag_candidate = range:: range_from_begin_end( w_gwas_begin, w_gwas_end )
                        ; ! tag_candidate.empty()
                        ;   tag_candidate.advance()
                        ) {
                    auto crps = tag_candidate.current_it().get_chrpos();
                    // Find the ref panel entries at the exact same chr:pos
                    auto ref_candidates = range:: range_from_begin_end(
                             std:: lower_bound( all_nearby_ref_data.begin() , all_nearby_ref_data.end(), chrpos(crps) )
                            ,std:: upper_bound( all_nearby_ref_data.begin() , all_nearby_ref_data.end(), chrpos(crps) )
                            );
                    vector<double> one_tag_zs;
                    vector<RefRecord const *> one_tag_its;
                    for(; ! ref_candidates.empty(); ref_candidates.advance()) {
                        auto & current_ref = ref_candidates.front_ref();
                        assert( current_ref.pos == crps.pos );

                        auto dir = decide_on_a_direction( current_ref
                                             , tag_candidate .current_it() );
                        switch(dir) {
                            break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED:
                                one_tag_zs.push_back( -tag_candidate.current_it().get_z() );
                                one_tag_its.push_back( &current_ref        );
                            break; case which_direction_t:: DIRECTION_AS_IS             :
                                one_tag_zs.push_back(  tag_candidate.current_it().get_z() );
                                one_tag_its.push_back( &current_ref        );
                            break; case which_direction_t:: NO_ALLELE_MATCH             : ;
                        }
                    }
                    // Finally, we make the decision on whether to include this tag or not.
                    // A tag is useful if there is exactly one reference SNP that corresponds
                    // to it (where 'correspond' means that the positions and alleles match).
                    assert(one_tag_zs.size() == one_tag_its.size());
                    if(1==one_tag_zs.size()) {
                        for(auto z : one_tag_zs)
                            tag_zs_.push_back(z);
                        for(auto it : one_tag_its)
                            tag_its_.push_back(it);
                    }
                }
            }

            assert(tag_zs_.size() == tag_its_.size());

            int const number_of_tags = tag_zs_.size();
            if(number_of_tags == 0)
                continue; // to the next window

            // Now, find suitable targets - i.e. anything in the reference panel in the narrow window
            // But some SNPs will have to be dropped, depending on --impute.range and --impute.snps
            vector<RefRecord const*>                unk2_its;
            {
                auto w_ref_narrow_begin = std:: lower_bound (   all_nearby_ref_data.begin() ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_start});
                auto w_ref_narrow_end   = std:: lower_bound (   all_nearby_ref_data.begin() ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_end});
                for(auto it = w_ref_narrow_begin; it<w_ref_narrow_end; ++it) {

                    bool skip_this_target = skipper_target->skip_me(chrm, &*it);
                    if(skip_this_target)
                        continue;

                    auto const & z12_for_this_SNP = it->z12;
                    auto z12_minmax = minmax_element(z12_for_this_SNP.begin(), z12_for_this_SNP.end());
                    if  (*z12_minmax.first == *z12_minmax.second){
                        continue;   // no variation in the allele counts for this SNP within
                                    //the ref panel, therefore useless for imputation
                                    //(Actually, I think this might never happen, as the DISCARD rule will
                                    // take this into account. Anyway, no harm to leave this in.
                    }
                    unk2_its    .push_back( &*it       );
                }
            }


            // Next few lines are kind of boring, just printing some statistics.

            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;
            int  number_of_all_targets                          = unk2_its.size();

            if(number_of_all_targets == 0)
                continue;


            // We now have at least one tag, and at least one target, so we can proceed

            // TODO: I still am including all the tags among the targets - change this?

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << number_of_tags                                 << " # SNPs in both (i.e. useful as tags)\n";
            cout << setw(8) << number_of_all_targets                          << " # target SNPs (anything in narrow window, will include some tags)\n";

            // Next, actually look up all the relevant SNPs within the reference panel
            vector<vector<int>> genotypes_for_the_tags  = from:: vector(tag_its_) |view::map| [](RefRecord const *rrp) { return rrp->z12; } |action::collect;
            vector<vector<int>> genotypes_for_the_unks  = from:: vector(unk2_its) |view::map| [](RefRecord const *rrp) { return rrp->z12; } |action::collect;

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
            mvn:: VecCol        C_inv_zs    = solve_a_matrix (C, mvn:: make_VecCol(tag_zs_));
            mvn:: Matrix        c           = make_c_unkn_tags_matrix( genotypes_for_the_tags
                                                         , genotypes_for_the_unks
                                                         , tag_its_
                                                         , unk2_its
                                                         , options:: opt_lambda
                                                         );

            // Next two lines are the imputation, and its quality.
            auto c_Cinv_zs = mvn:: multiply_matrix_by_colvec_giving_colvec(c, C_inv_zs);
            auto   Cinv_c  =     mvn:: multiply_NoTrans_Trans(invert_a_matrix (C) , c);

            assert( (int)Cinv_c.size1() == number_of_tags );
            assert( (int)Cinv_c.size2() == number_of_all_targets );


            // Finally, print out the imputations
            assert(number_of_all_targets == ssize(unk2_its));
            mvn:: Matrix to_store_one_imputation_quality(1,1); // a one-by-one-matrix
            for(int i=0; i<number_of_all_targets; ++i) {
                    auto && target = unk2_its.at(i);
                    auto pos  = target->pos;
                    auto SNPname = target->ID;

                    auto imp_qual = [&](){ // multiply one vector by another, to directly get
                                            // the diagonal of the imputation quality matrix.
                                            // This should be quicker than fully computing
                                            // c' * inv(C) * c
                        auto lhs = gsl_matrix_const_submatrix(     c.get(), i, 0, 1, number_of_tags);
                        auto rhs = gsl_matrix_const_submatrix(Cinv_c.get(), 0, i, number_of_tags, 1);
                        // TODO: Maybe move the next line into mvn.{hh,cc}?
                        const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &lhs.matrix, &rhs.matrix, 0, to_store_one_imputation_quality.get());
                        assert(res_0 == 0);
                        return to_store_one_imputation_quality(0,0);
                    }();

                    (*out_stream_ptr)
                        << "chr" << chrm << ':' << pos
                        << '\t' << c_Cinv_zs(i)
                        << '\t' << SNPname
                        << '\t' << target->ref
                        << '\t' << target->alt
                        << '\t' << imp_qual
                        << endl;
            }
        }
    }
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
        , vector<RefRecord const *              > const & tag_its
        , vector<RefRecord const *              > const & unk_its
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
            ,   RefRecord const *
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

                if(tag_its_k == unk_its_u) { // identical pointer value, i.e. same reference entry
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
    ( RefRecord                       const & r
    , SNPiterator     const & g
    ) {
        auto const rp_ref = r.ref;
        auto const rp_alt = r.alt;
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
