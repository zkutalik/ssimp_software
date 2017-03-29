#include <iostream>
#include <limits>
#include <cassert>

#include "options.hh"
#include "file.reading.hh"

#include "other/DIE.hh"
#include "other/PP.hh"

using std:: cout;
using std:: endl;

#include "fwd/src/ssimp.hh"

int main(int argc, char **argv) {
    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);

    cout.imbue(std::locale("")); // apply the user's locale, for example the thousands separator

    if(!options:: opt_raw_ref.empty()) {
        PP(options:: opt_raw_ref);
        auto raw_ref_file = file_reading:: read_in_a_raw_ref_file(options:: opt_raw_ref);
        PP(raw_ref_file->number_of_snps());

        ssimp:: quickly_list_the_regions(raw_ref_file);
    }
}

namespace ssimp{

FWD(ssimp)
void quickly_list_the_regions(file_reading:: GenotypeFileHandle raw_ref_file) {
    int total_number_of_SNPs = raw_ref_file->number_of_snps();
    PP (total_number_of_SNPs);

    auto const b = file_reading:: SNPiterator:: begin_from_file(raw_ref_file);
    auto const e = file_reading:: SNPiterator::   end_from_file(raw_ref_file);

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

} // namespace ssimp
