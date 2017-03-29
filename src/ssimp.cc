#include <iostream>

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

    auto b = file_reading:: SNPiterator:: begin_from_file(raw_ref_file);
    auto e = file_reading:: SNPiterator::   end_from_file(raw_ref_file);

    while(b <  e) {
        PP(b.get_chrpos());
        ++b;
    }
    PP("end");
}

} // namespace ssimp
