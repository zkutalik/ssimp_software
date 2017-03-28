#include <iostream>

#include "options.hh"

#include "other/DIE.hh"
#include "other/PP.hh"

using std:: cout;
using std:: endl;

int main(int argc, char **argv) {
    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);

    PP(options:: opt_raw_ref);
}
