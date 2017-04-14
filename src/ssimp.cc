#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>

#include "options.hh"
#include "file.reading.hh"

#include "other/utils.hh"
#include "other/mvn.hh"
#include "other/DIE.hh"
#include "other/PP.hh"

using std:: cout;
using std:: endl;
using std:: string;
using std:: vector;
using std:: setw;

using file_reading:: chrpos;
using file_reading:: SNPiterator;
using file_reading:: GenotypeFileHandle;
using file_reading:: GwasFileHandle;

using utils:: ssize;

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
vector<vector<int>>
lookup_genotypes( vector<chrpos>                     const &  SNPs_in_the_intersection
                , file_reading:: CacheOfRefPanelData       &  cache
                );
static
mvn:: SquareMatrix
make_C_tag_tag_matrix_and_invert(
                    vector<vector<int>>              const & genotypes_for_the_tags );
static
mvn:: VecCol
make_C_tag_tag_matrix_and_solve(
                    vector<vector<int>>              const & genotypes_for_the_tags
                  , mvn:: VecCol                             zs
                    );
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int>>              const & genotypes_for_the_tags
        , vector<vector<int>>              const & genotypes_for_the_unks
        );

enum class which_direction_t { DIRECTION_SHOULD_BE_REVERSED
                             , NO_ALLELE_MATCH
                             , DIRECTION_AS_IS };
static
    which_direction_t decide_on_a_direction
    ( SNPiterator<GenotypeFileHandle> const & //r
    , SNPiterator<GwasFileHandle>     const & //g
    );
} // namespace ssimp

int main(int argc, char **argv) {
    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);

    cout.imbue(std::locale("")); // apply the user's locale, for example the thousands separator

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
void impute_all_the_regions( file_reading:: GenotypeFileHandle         ref_panel
                             , file_reading:: GwasFileHandle             gwas
                             ) {
    auto const b_ref  = begin_from_file(ref_panel);
    auto const e_ref  =   end_from_file(ref_panel);
    auto const b_gwas = begin_from_file(gwas);
    auto const e_gwas =   end_from_file(gwas);

    PP(options:: opt_window_width
      ,options:: opt_flanking_width
            );

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
            auto w_ref_narrow_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start});
            auto w_ref_narrow_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end  });
            if(w_ref_narrow_begin == c_end)
                break; // Finished with this chromosome
            if(w_ref_narrow_begin == w_ref_narrow_end)
                continue; // Nothing to impute, skip to the next region

            auto w_ref_wide_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_ref_wide_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end   + options:: opt_flanking_width });

            // Look up this region in the GWAS (taking account of the flanking width also)
            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            // Which SNPs are in both *wide* windows, i.e. useful as tag SNPs
            vector<chrpos> SNPs_in_the_intersection;
            vector<double> zs_for_the_tags;
            {
                auto r = w_ref_wide_begin;
                auto g = w_gwas_begin;
                while(r < w_ref_wide_end && g < w_gwas_end) {
                    if(r.get_chrpos() == g.get_chrpos()) {
                        { // one last check, that the alleles match up
                            which_direction_t which_dir = decide_on_a_direction(r, g);
                            switch(which_dir) {
                                break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED:
                                    {
                                        auto z_in_ref_direction = -g.get_z();
                                        (void)z_in_ref_direction;
                                        SNPs_in_the_intersection.push_back(r.get_chrpos());
                                        zs_for_the_tags.push_back(z_in_ref_direction);
                                    }
                                break; case which_direction_t:: DIRECTION_AS_IS:
                                    {
                                        auto z_in_ref_direction =  g.get_z();
                                        (void)z_in_ref_direction;
                                        SNPs_in_the_intersection.push_back(r.get_chrpos());
                                        zs_for_the_tags.push_back(z_in_ref_direction);
                                    }
                                break; case which_direction_t:: NO_ALLELE_MATCH: ; // Just ignore this
                            }
                        }
                        ++g; // because g might have the same position multiple times, but the reference panel shouldn't
                             // TODO: check that the reference panel doesn't indeed have unique chrpos.
                    }
                    else if(r.get_chrpos() <  g.get_chrpos()) {
                        ++r; // skip this one
                        r = std:: lower_bound(r, w_ref_wide_end, g.get_chrpos());
                    }
                    else if(r.get_chrpos() >  g.get_chrpos()) {
                        ++g; // skip this one
                        g = std:: lower_bound(g, w_gwas_end, r.get_chrpos());
                    }
                }
            }
            assert(SNPs_in_the_intersection.size() == zs_for_the_tags.size());

            int const number_of_tags = SNPs_in_the_intersection.size();
            if(number_of_tags == 0)
                continue;
            assert(number_of_tags > 0);

            // Now, find suitable targets - i.e. anything in the reference panel in the narrow window
            vector<chrpos>  SNPs_all_targets;
            for(auto it = w_ref_narrow_begin; it<w_ref_narrow_end; ++it) {
                // actually, we should think about ignoring SNPs in certain situations
                auto allele_alt =it.get_allele_alt();
                auto has_more_than_one_alt_allele = allele_alt.find(',') != std::string::npos;
                if(has_more_than_one_alt_allele)
                    continue;

                auto const & z12_for_this_SNP = cache.lookup_one_chr_pos(it.get_chrpos());
                auto z12_minmax = minmax_element(z12_for_this_SNP.begin(), z12_for_this_SNP.end());
                if  (*z12_minmax.first == *z12_minmax.second){
                    continue; // no variation in this SNP within the ref panel, therefore useless for imputation
                }
                SNPs_all_targets.push_back( it.get_chrpos() );
            }


            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;
            int  number_of_all_targets                          = SNPs_all_targets.size();

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << number_of_tags                                 << " # SNPs in both (i.e. useful as tags)\n";
            cout << setw(8) << number_of_all_targets                          << " # target SNPs (anything in narrow window, will include some tags)\n";
            auto genotypes_for_the_tags = lookup_genotypes( SNPs_in_the_intersection, cache );
            auto genotypes_for_the_unks = lookup_genotypes( SNPs_all_targets        , cache );
            assert(number_of_tags == utils:: ssize(genotypes_for_the_tags));
            assert(number_of_all_targets == utils:: ssize(genotypes_for_the_unks));

            static_assert( std:: is_same< vector<vector<int>> , decltype(genotypes_for_the_tags) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation
            static_assert( std:: is_same< vector<vector<int>> , decltype(genotypes_for_the_unks) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation

            int const N_ref = genotypes_for_the_tags.at(0).size();
            assert(N_ref > 0);

            mvn:: SquareMatrix C_inv = make_C_tag_tag_matrix_and_invert(genotypes_for_the_tags);
            mvn:: VecCol C_inv_zs = make_C_tag_tag_matrix_and_solve(genotypes_for_the_tags, mvn:: make_VecCol(zs_for_the_tags));
            auto old_way =
                    multiply_matrix_by_colvec_giving_colvec
                                              ( C_inv
                                              , mvn:: make_VecCol(zs_for_the_tags)
                                              );
            {
                auto diff = old_way - C_inv_zs;
                double mn, mx;
                gsl_vector_minmax(diff.get(), &mn, &mx);
                assert(mn > -1e-10);
                assert(mx <  1e-10);
            }
            mvn:: Matrix      c = make_c_unkn_tags_matrix(
                    genotypes_for_the_tags
                    ,genotypes_for_the_unks
                    );

            auto cCzs = multiply_matrix_by_colvec_giving_colvec
                        ( c
                        , multiply_matrix_by_colvec_giving_colvec
                          ( C_inv
                          , mvn:: make_VecCol(zs_for_the_tags)
                          )
                        );
            PP(cCzs);
            assert(number_of_all_targets == ssize(cCzs));
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
        auto rel = m.insert( std:: make_pair(b.get_SNPname(), b.get_chrpos()) );
        rel.second || DIE("same SNPname twice in the ref panel [" << b.get_SNPname() << "]");
    }
    return m;
}
static
vector<vector<int>>
lookup_genotypes( vector<chrpos>                     const &  snps
                , file_reading:: CacheOfRefPanelData       &  cache
                ) {
    vector<vector<int>> many_calls;
    for(auto crps : snps) {
        many_calls.push_back( cache.lookup_one_chr_pos(crps));
    }
    assert(many_calls.size() == snps.size());
    return many_calls;
}

static
mvn:: SquareMatrix
make_C_tag_tag_matrix_and_invert(
                    vector<vector<int>>              const & genotypes_for_the_tags ) {
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
                        if(c_kl > 1.0) {
                            assert(c_kl-1.0 < 1e-5);
                            c_kl = 1.0;
                        }
                        assert(c_kl >= -1.0);
                        assert(c_kl <=  1.0);
                        if(k==l)
                            assert(c_kl == 1.0);
                        else
                            assert(c_kl <  1.0);
                        C.set(k,l,c_kl);
                    }
                }
                return invert_a_matrix(std::move(C));
}
static
mvn:: VecCol
make_C_tag_tag_matrix_and_solve(
                    vector<vector<int>>              const & genotypes_for_the_tags
                  , mvn:: VecCol                             zs
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
                        if(c_kl > 1.0) {
                            assert(c_kl-1.0 < 1e-5);
                            c_kl = 1.0;
                        }
                        assert(c_kl >= -1.0);
                        assert(c_kl <=  1.0);
                        if(k==l)
                            assert(c_kl == 1.0);
                        else
                            assert(c_kl <  1.0);
                        C.set(k,l,c_kl);
                    }
                }
                return solve_a_matrix(std::move(C), zs);
}
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int>>              const & genotypes_for_the_tags
        , vector<vector<int>>              const & genotypes_for_the_unks
        ) {
    int const number_of_tags = genotypes_for_the_tags.size();
    int const number_of_all_targets = genotypes_for_the_unks.size();
    int const N_ref = genotypes_for_the_tags.at(0).size();
    assert(N_ref > 0);

    mvn:: Matrix      c(number_of_all_targets, number_of_tags);
    for(int k=0; k<number_of_tags; ++k) {
        for(int u=0; u<number_of_all_targets; ++u) {
            assert(N_ref        == utils:: ssize(genotypes_for_the_tags.at(k)));
            assert(N_ref        == utils:: ssize(genotypes_for_the_unks.at(u)));
            double c_ku = gsl_stats_int_correlation( &genotypes_for_the_tags.at(k).front(), 1
                                                   , &genotypes_for_the_unks.at(u).front(), 1
                                                   , N_ref );
            if(c_ku > 1.0) {
                assert(c_ku-1.0 < 1e-5);
                c_ku = 1.0;
            }
            assert(c_ku >= -1.0);
            assert(c_ku <=  1.0);
            c.set(u,k,c_ku);
        }
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
