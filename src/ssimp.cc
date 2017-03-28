#include <iostream>

#include "options.hh"

#include "other/DIE.hh"

using std:: cout;
using std:: endl;

int main(int argc, char **argv) {
    // all options now read. Start checking they are all present
    options:: read_in_all_command_line_options(argc, argv);
}
