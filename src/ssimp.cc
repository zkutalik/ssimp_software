#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_statistics_int.h>

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

namespace ssimp {
// Some forward declarations
static
void quickly_list_the_regions( file_reading:: GenotypeFileHandle         raw_ref_file
                             , file_reading:: GwasFileHandle             gwas
                             );
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
        ssimp:: quickly_list_the_regions(raw_ref_file, gwas);
        ssimp:: impute_all_the_regions(raw_ref_file, gwas);
    }
}

namespace ssimp{

using file_reading:: GenotypeFileHandle;
using file_reading:: GwasFileHandle;

void quickly_list_the_regions( file_reading:: GenotypeFileHandle         ref_panel
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

        // First, find the begin and end of this chromosome
        auto c_begin = std:: lower_bound(b_ref, e_ref, chrpos{chrm, std::numeric_limits<int>::lowest() });
        auto c_end   = std:: lower_bound(b_ref, e_ref, chrpos{chrm, std::numeric_limits<int>::max()  });
        assert(c_end >= c_begin);
        if(c_begin != c_end)
            assert(c_begin.get_chrpos().pos >= 0); // first position is at least zero

        for(int w = 0; ; ++w ) {
            int current_window_start = w     * options:: opt_window_width;
            int current_window_end   = (w+1) * options:: opt_window_width;
            auto w_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start});
            auto w_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end  });
            if(w_begin == c_end)
                break; // Finished with this chromosome
            if(w_begin == w_end)
                continue; // Empty region, just skip it

            // Look up this region in the GWAS (taking account of the flanking width also)
            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            // Which SNPs are in both, i.e. useful as tag SNPs
            vector<chrpos> SNPs_in_the_intersection;
            {
                auto r = w_begin;
                auto g = w_gwas_begin;
                while(r < w_end && g < w_gwas_end) {
                    if(r.get_chrpos() == g.get_chrpos()) {
                        SNPs_in_the_intersection.push_back(r.get_chrpos());
                        ++g; // because g might have the same position multiple times, but the reference panel shouldn't
                             // TODO: check that the reference panel doesn't indeed have unique chrpos.
                    }
                    else if(r.get_chrpos() <  g.get_chrpos()) {
                        ++r; // skip this one
                        r = std:: lower_bound(r, w_end, g.get_chrpos());
                    }
                    else if(r.get_chrpos() >  g.get_chrpos()) {
                        ++g; // skip this one
                        g = std:: lower_bound(g, w_gwas_end, r.get_chrpos());
                    }
                }
            }
            { // verify. Feel free to just delete these coming lines
                vector<chrpos> SNPs_in_the_intersection_;
                std:: set_intersection( w_begin     , w_end
                                      , w_gwas_begin, w_gwas_end
                                      , std:: back_inserter( SNPs_in_the_intersection_ ));
                assert(SNPs_in_the_intersection == SNPs_in_the_intersection_);
            }


            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_ref_panel_in_this_region = w_end      - w_begin;
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_ref_panel_in_this_region << " # RefPanel SNPs in this window\n";
            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << SNPs_in_the_intersection.size()                << " # SNPs in both (i.e. useful as tags)\n";
        }
    }

}
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
            auto w_begin = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_start});
            auto w_end   = std:: lower_bound(c_begin, c_end, chrpos{chrm,current_window_end  });
            if(w_begin == c_end)
                break; // Finished with this chromosome
            if(w_begin == w_end)
                continue; // Empty region, just skip it

            // Look up this region in the GWAS (taking account of the flanking width also)
            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            // Which SNPs are in both, i.e. useful as tag SNPs
            vector<chrpos> SNPs_in_the_intersection;
            {
                auto r = w_begin;
                auto g = w_gwas_begin;
                while(r < w_end && g < w_gwas_end) {
                    if(r.get_chrpos() == g.get_chrpos()) {
                        SNPs_in_the_intersection.push_back(r.get_chrpos());
                        ++g; // because g might have the same position multiple times, but the reference panel shouldn't
                             // TODO: check that the reference panel doesn't indeed have unique chrpos.
                    }
                    else if(r.get_chrpos() <  g.get_chrpos()) {
                        ++r; // skip this one
                        r = std:: lower_bound(r, w_end, g.get_chrpos());
                    }
                    else if(r.get_chrpos() >  g.get_chrpos()) {
                        ++g; // skip this one
                        g = std:: lower_bound(g, w_gwas_end, r.get_chrpos());
                    }
                }
            }
            { // verify. Feel free to just delete these coming lines
                vector<chrpos> SNPs_in_the_intersection_;
                std:: set_intersection( w_begin     , w_end
                                      , w_gwas_begin, w_gwas_end
                                      , std:: back_inserter( SNPs_in_the_intersection_ ));
                assert(SNPs_in_the_intersection == SNPs_in_the_intersection_);
            }

            int const number_of_tags = SNPs_in_the_intersection.size();
            if(number_of_tags == 0)
                continue;
            assert(number_of_tags > 0);


            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_ref_panel_in_this_region = w_end      - w_begin;
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_ref_panel_in_this_region << " # RefPanel SNPs in this window\n";
            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << number_of_tags                                 << " # SNPs in both (i.e. useful as tags)\n";
            auto genotypes_for_the_tags = lookup_genotypes( SNPs_in_the_intersection, cache );

            static_assert( std:: is_same< vector<vector<int>> , decltype(genotypes_for_the_tags) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation

            int const N_ref = genotypes_for_the_tags.at(0).size();
            assert(N_ref > 0);
            PP(N_ref);


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
                    PP(k,l,c_kl);
                    assert(c_kl >= -1.0);
                    assert(c_kl <=  1.0);
                    if(k==l)
                        assert(c_kl == 1.0);
                    else
                        assert(c_kl <  1.0);
                }
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

} // namespace ssimp
