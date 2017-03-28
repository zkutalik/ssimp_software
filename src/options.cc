#include "options.hh"

#include <getopt.h>
#include <iostream>

#include "other/DIE.hh"

namespace options {

void read_in_all_command_line_options(int argc, char **argv) {
    while(1) { // while there are still more options to be processed
        int long_option_index;
        static struct option long_options[] = {
            //{"header.in",  no_argument, &opt_headerin,  1 }, // a boolean flag
            {"raw.ref",  required_argument, 0,  2 }, // last arg on each line must be greater than 1
            {0,         0,                 0,  0 } // must have this line of zeroes at the end
        };
        int c = getopt_long(argc, argv, "-", long_options, &long_option_index);
        if (c == -1)
            break;
        if (c == 1) { // non-option
            DIE("There shouldn't be any non-options args to this program: ["<<optarg<<"]");
        }
        if (c == 2) {
            //assert(string("z_compare") == long_options[long_option_index].name);
            //zs_to_compare.emplace_back( optarg, read_in_a_vector(optarg) );
        }
    }
}

} // namespace options
