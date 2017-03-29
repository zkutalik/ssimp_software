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

    if(!options:: opt_raw_ref.empty()) {
        PP(options:: opt_raw_ref);
        auto raw_ref_file = file_reading:: read_in_a_raw_ref_file(options:: opt_raw_ref);
        PP(raw_ref_file->number_of_snps());
    }
}
