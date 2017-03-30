#include <iostream>
#include <limits>
#include <cassert>
#include <unordered_map>

#include "options.hh"
#include "file.reading.hh"

#include "other/DIE.hh"
#include "other/PP.hh"

using std:: cout;
using std:: endl;
using std:: string;

namespace ssimp {
// Some forward declarations
void quickly_list_the_regions(file_reading:: GenotypeFileHandle raw_ref_file);
std:: unordered_map<string, file_reading:: chrpos>
            map_rs_to_chrpos( file_reading:: GenotypeFileHandle raw_ref_file);
} // namespace ssimp

int main(int argc, char **argv) {
    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);

    cout.imbue(std::locale("")); // apply the user's locale, for example the thousands separator

    if(!options:: opt_raw_ref.empty() && !options:: opt_gwas_filename.empty()) {
        PP( options:: opt_raw_ref
          , options:: opt_gwas_filename);

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

        PP(  gwas->number_of_snps());
        PP(raw_ref_file->number_of_snps());

        ssimp:: quickly_list_the_regions(raw_ref_file);
    }
}

namespace ssimp{

using file_reading:: GenotypeFileHandle;
using file_reading:: GwasFileHandle;

void quickly_list_the_regions(file_reading:: GenotypeFileHandle raw_ref_file) {
    int total_number_of_SNPs = raw_ref_file->number_of_snps();
    PP (total_number_of_SNPs);

    auto const b = file_reading:: SNPiterator<GenotypeFileHandle>:: begin_from_file(raw_ref_file);
    auto const e = file_reading:: SNPiterator<GenotypeFileHandle>::   end_from_file(raw_ref_file);

    /*
     * 16050075
     * 16050115
     * 16050213
     * ....
     * 16644605
     * 16644612
     * 16644620
     */

    PP(options:: opt_window_width);

    for(int chrm = 22; chrm <= 22; ++chrm) {

        // First, find the begin and end of this chromosome
        auto c_begin = std:: lower_bound(b, e, file_reading:: chrpos{chrm, 0 });
        auto c_end   = std:: lower_bound(b, e, file_reading:: chrpos{chrm, std::numeric_limits<int>::max()  });
        assert(c_end >= c_begin);

        for(int w = 0; ; ++w ) {
            int current_window_start = w     * options:: opt_window_width;
            int current_window_end   = (w+1) * options:: opt_window_width;
            auto w_begin = std:: lower_bound(c_begin, c_end, file_reading:: chrpos{chrm,current_window_start});
            auto w_end   = std:: lower_bound(w_begin, c_end, file_reading:: chrpos{chrm,current_window_end  });
            if(w_begin == c_end)
                break; // Finished with this chromosome
            if(w_begin == w_end)
                continue; // Empty region, just skip it
            PP(chrm, w, current_window_start, current_window_end-1, w_end - w_begin);
            // We have at least one SNP here, let's print some stuff about it
            PP( *w_begin
                    , w_begin.get_SNPname()
                    , w_begin.get_allele_ref()
                    , w_begin.get_allele_alt()
                    );
        }
    }

}

std:: unordered_map<string, file_reading:: chrpos>
            map_rs_to_chrpos( file_reading:: GenotypeFileHandle raw_ref_file ) {
    auto       b = file_reading:: SNPiterator<GenotypeFileHandle>:: begin_from_file(raw_ref_file);
    auto const e = file_reading:: SNPiterator<GenotypeFileHandle>::   end_from_file(raw_ref_file);
    std:: unordered_map<string, file_reading:: chrpos> m;
    for(;b<e; ++b) {
        auto rel = m.insert( std:: make_pair(b.get_SNPname(), b.get_chrpos()) );
        rel.second || DIE("same SNPname twice in the ref panel [" << b.get_SNPname() << "]");
    }
    PP(m.size());
    return m;
}

} // namespace ssimp
