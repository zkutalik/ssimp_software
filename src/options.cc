#include "options.hh"

#include <unistd.h> // for 'unlink'
#include <getopt.h>
#include <iostream>
#include <fstream>

#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "bits.and.pieces/utils.hh"

#include "file.reading.hh"
#include "range/range_view.hh"
#include "range/range_action.hh"
#include "format/format.hh"

using std:: string;

namespace view   = range :: view;
namespace action = range :: action;

namespace options {

        std::vector<std:: string>           opt_non_options; // if there are two or three of these, copy them into opt_gwas_filename, opt_out, and opt_raw_ref
        std:: string            opt_raw_ref;
        std:: string            opt_gwas_filename;
        std:: string            opt_out;
        std:: string            opt_log; // copy of whatever is sent to the console.
        int                     opt_window_width = 1'000'000;
        int                     opt_flanking_width = 250'000;
        std:: string            opt_lambda  = "2/sqrt(n)";

        std:: string            opt_impute_range;
        std:: string            opt_impute_snps;
        std::unique_ptr<std::unordered_set<std::string>>    opt_impute_snps_as_a_uset;
        double                  opt_impute_maf =0.0;

        // The next few are like the --impute.* above, but applying to tags instead
        std:: string            opt_tag_range;
        std:: string            opt_tag_snps;
        std::unique_ptr<std::unordered_set<std::string>>    opt_tag_snps_as_a_uset;
        double                  opt_tag_maf =0.0; // target not imputed unless maf (in reference) is at least this.

        bool                    opt_reimpute_tags = false;
        std:: string            opt_tags_used_output;

        std:: string            opt_sample_names;
        opt_missingness_t       opt_missingness;


        std:: vector<std::function<void(void)>>    list_of_tasks_to_run_at_exit;

void read_in_all_command_line_options(int argc, char **argv) {
    while(1) { // while there are still more options to be processed
        int long_option_index;
        static struct option long_options[] = {
            //{"header.in",  no_argument, &opt_headerin,  1 }, // a boolean flag
            {"ref"                ,  required_argument, 0,  2 }, // last arg on each line must be greater than 1
            {"window.width"       ,  required_argument, 0,  3 },
            {"flanking.width"     ,  required_argument, 0,  4 },
            {"gwas"               ,  required_argument, 0,  5 },
            {"lambda"             ,  required_argument, 0,  6 },
            {"out"                ,  required_argument, 0,  7 },
            {"impute.range"       ,  required_argument, 0,  8 },
            {"impute.snps"        ,  required_argument, 0,  9 },
            {"impute.maf"         ,  required_argument, 0, 10 },
            {"tag.range"          ,  required_argument, 0, 11 },
            {"tag.snps"           ,  required_argument, 0, 12 },
            {"tag.maf"            ,  required_argument, 0, 13 },
            {"reimpute.tags"      ,        no_argument, 0, 14 }, // one-by-one, reimpute each tag by masking it
            {"sample.names"       ,  required_argument, 0, 15 },
            {"tags.used.output"   ,  required_argument, 0, 16 },
            {"log"                ,  required_argument, 0, 17 },
            {"missingness"        ,  required_argument, 0, 18 },
            {0                    ,  0                , 0,  0 } // must have this line of zeroes at the end
        };
        int c = getopt_long(argc, argv, "-", long_options, &long_option_index);
        if (c == '?')
            DIE("problem with options");
        if (c == -1)
            break;
        if (c == 1) { // non-option
            opt_non_options.push_back(optarg);
            // Save all the non-option arguments. There should be two or three of them
        }
        if (c == 2) {
            assert(string("ref") == long_options[long_option_index].name);
            options:: opt_raw_ref = optarg;
        }
        if (c == 3) {
            options::  opt_window_width  == 1'000'000 || DIE("--window width specified twice?");
            assert(string("window.width") == long_options[long_option_index].name);
            options::  opt_window_width  = utils:: lexical_cast<int>(optarg);
        }
        if (c == 4) {
            assert(string("flanking.width") == long_options[long_option_index].name);
            options::  opt_flanking_width  = utils:: lexical_cast<int>(optarg);
        }
        if (c == 5) {
            assert(string("gwas") == long_options[long_option_index].name);
            options::  opt_gwas_filename = optarg;
        }
        if (c == 6) {
            assert(string("lambda") == long_options[long_option_index].name);
            options::  opt_lambda  = optarg;
        }
        if (c == 7) {
            assert(string("out") == long_options[long_option_index].name);
            options::  opt_out = optarg;
        }
        if (c == 8) {
            options:: opt_impute_range.empty() || DIE("--impute.range specified twice?");
            assert(string("impute.range") == long_options[long_option_index].name);
            options::  opt_impute_range = optarg;
        }
        if (c == 9) {
            options:: opt_impute_snps.empty() || DIE("--impute.snps specified twice?");
            assert(string("impute.snps") == long_options[long_option_index].name);
            options::  opt_impute_snps = optarg;
        }
        if (c == 10) {
            options:: opt_impute_maf == 0.0 || DIE("--impute.maf specified twice?");
            assert(string("impute.maf") == long_options[long_option_index].name);
            options::  opt_impute_maf = utils:: lexical_cast<double>(optarg);
        }
        if (c == 11) {
            options:: opt_tag_range.empty() || DIE("--tag.range specified twice?");
            assert(string("tag.range") == long_options[long_option_index].name);
            options::  opt_tag_range = optarg;
        }
        if (c == 12) {
            options:: opt_tag_snps.empty() || DIE("--tag.snps specified twice?");
            assert(string("tag.snps") == long_options[long_option_index].name);
            options::  opt_tag_snps = optarg;
        }
        if (c == 13) {
            options:: opt_tag_maf == 0.0 || DIE("--tag.maf specified twice?");
            assert(string("tag.maf") == long_options[long_option_index].name);
            options::  opt_tag_maf = utils:: lexical_cast<double>(optarg);
        }
        if (c == 14) {
            assert(string("reimpute.tags") == long_options[long_option_index].name);
            options::  opt_reimpute_tags = true;
        }
        if (c == 15) {
            options:: opt_sample_names.empty() || DIE("--sample.names specified twice?");
            assert(string("sample.names") == long_options[long_option_index].name);
            options::  opt_sample_names = optarg;
        }
        if (c == 16) {
            options:: opt_tags_used_output.empty() || DIE("--tags.used.output specified twice?");
            assert(string("tags.used.output") == long_options[long_option_index].name);
            options::  opt_tags_used_output = optarg;
        }
        if (c == 17) {
            options:: opt_log.empty() || DIE("--log specified twice?");
            assert(string("log") == long_options[long_option_index].name);
            options::  opt_log = optarg;
        }
        if (c == 18) {
            assert(string("missingness") == long_options[long_option_index].name);
            if(false) {}
            else if(string("none") == optarg) {
                options::  opt_missingness = opt_missingness_t:: NAIVE;
            }
            else if(string("dep") == optarg) {
                options::  opt_missingness = opt_missingness_t:: DEPENDENCY_MAXIMUM;
            }
            else if(string("ind") == optarg) {
                options::  opt_missingness = opt_missingness_t:: INDEPENDENCE;
            }
            else {
                DIE("--missingness argument ([" << optarg << "]) not understood. Valid values: dep,ind,none");
            }
        }
    }
}

void adjust_sample_names_if_it_is_magical() { // check the filename exists. If not, treat it as a magical filename
    std:: ifstream test_if_file_exists(options::  opt_sample_names.c_str());
    if(!test_if_file_exists) {
        // maybe everything after the second-last '/' is a filter?
        size_t lastslash = options::  opt_sample_names.find_last_of('/');
        if(lastslash != string::npos) {
            size_t secondlastslash = options::  opt_sample_names.find_last_of('/', lastslash-1);
            if(secondlastslash != string:: npos) {
                string underlying_filename = options::  opt_sample_names.substr(0,secondlastslash);
                string sample_field        = options::  opt_sample_names.substr(secondlastslash+1, lastslash-secondlastslash-1);
                string filter              = options::  opt_sample_names.substr(lastslash+1);
                auto   filter_field_pair   = utils:: tokenize(filter, '=');
                filter_field_pair.size() == 2 || DIE("Should be just one '=' in [" << filter << "]");
                string filter_field        = filter_field_pair.at(0);
                string filter_value        = filter_field_pair.at(1);
                std:: ifstream full_panel (underlying_filename);
                full_panel || DIE("Can't find the full panel file [" << underlying_filename << "]");

                string header;
                getline(full_panel, header);
                full_panel || DIE("Couldn't read in a header line from [" << underlying_filename << "]");

                char delimiter = file_reading:: decide_delimiter(header);
                auto header_fields = utils:: tokenize( header, delimiter);

                int offset_of_sample_field = -1;
                int offset_of_filter_field = -1;

                using utils:: operator<<;

                view:: enumerate_vector(header_fields)
                |action:: unzip_foreach|
                [&](int i, string const & s) {
                    if( s == sample_field) {
                        offset_of_sample_field == -1 || DIE("double field name in [" << header_fields << "]?");
                        offset_of_sample_field = i;
                    }
                    if( s == filter_field) {
                        offset_of_filter_field == -1 || DIE("double field name in [" << header_fields << "]?");
                        offset_of_filter_field = i;
                    }
                };

                offset_of_sample_field != -1 || DIE("Couldn't find ["<<sample_field<<"] field in [" << header_fields << "]?");
                offset_of_filter_field != -1 || DIE("Couldn't find ["<<filter_field<<"] field in [" << header_fields << "]?");

                std:: vector<string> filtered_samples_to_use;
                string dataline;
                while(getline(full_panel, dataline)) {
                    auto dataline_split = utils:: tokenize(dataline, delimiter);
                    if(dataline_split.at(offset_of_filter_field) != filter_value)
                        continue;
                    filtered_samples_to_use.push_back(dataline_split.at(offset_of_sample_field));
                }

                char temporary_dirname[] = "/tmp/ssimp_XXXXXX";
                mkdtemp(temporary_dirname) || DIE("Couldn't create temporary filename for use with --sample_names [" << temporary_dirname << "]" );
                auto temporary_filename = AMD_FORMATTED_STRING("{0}/sample.names", temporary_dirname);
                std:: ofstream file_of_filtered_samples_to_use(temporary_filename);
                file_of_filtered_samples_to_use || DIE("Couldn't create temporary filename for use with --sample_names");

                options:: list_of_tasks_to_run_at_exit.push_back( // ensure this little task is executed at program exit
                    [temporary_dirname, temporary_filename](void) -> void{
                        int ret1 = unlink(temporary_filename.c_str());
                        ret1 == 0 || DIE("Couldn't delete temporary filename [" << temporary_filename << "]");
                        int ret2 =  rmdir(temporary_dirname);
                        ret2 == 0 || DIE("Couldn't delete temporary dirname ["  << temporary_dirname << "]");
                    });

                for(auto && one_sample : filtered_samples_to_use) {
                    file_of_filtered_samples_to_use << one_sample << '\n';
                }
                file_of_filtered_samples_to_use.close();
                options::  opt_sample_names = temporary_filename;
            }
        }
    }
}

} // namespace options
