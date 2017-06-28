#include <unistd.h> // for 'chdir' and 'unlink'
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

#include <gzstream.h>

#include "options.hh"
#include "file.reading.hh"

#include "usage_special.hh"

#include "bits.and.pieces/utils.hh"
#include "mvn/mvn.hh"
#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"
#include "range/range_view.hh"
#include "range/range_from.hh"
#include "range/range_action.hh"
#include "format/format.hh"

#include "file.reading.vcfgztbi.hh"

#include "logging.hh"

namespace view = range:: view;
namespace from = range:: from;
namespace action = range:: action;

using std:: cout;
using std:: endl;
using std:: string;
using std:: vector;
using std:: setw;
using std:: unique_ptr;
using std:: make_unique;
using std:: ofstream;
using std:: ostream;

using file_reading:: chrpos;
using file_reading:: SNPiterator;

using utils:: ssize;
using utils:: print_type;
using utils:: stdget0;
using utils:: stdget1;

using tbi:: RefRecord;


namespace ssimp {
// Some forward declarations

static
void impute_all_the_regions(   string                                   filename_of_vcf
                             , file_reading:: GwasFileHandle_NONCONST   gwas
                             );
static
mvn::SquareMatrix apply_lambda_square (mvn:: SquareMatrix copy, double lambda);
static
mvn:: Matrix apply_lambda_rect   (   mvn:: Matrix copy
                        ,   vector<RefRecord const *> const & unk2_its
                        ,   vector<RefRecord const *> const & tag_its_
                        , double lambda);

static
int compute_number_of_effective_tests_in_C_nolambda(mvn:: SquareMatrix const & C_smalllambda);
static
double  adjusted_imputation_quality (double iq_1e8_simple, int N_ref, int number_of_effective_tests);

static
mvn:: SquareMatrix
make_C_tag_tag_matrix( vector<vector<int> const *>      const & genotypes_for_the_tags);
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int> const *>      const & genotypes_for_the_tags
        , vector<vector<int> const *>      const & genotypes_for_the_unks
        , vector<RefRecord const *              > const & tag_its
        , vector<RefRecord const *              > const & unk_its
        );
using reimputed_tags_in_this_window_t = std:: unordered_map< RefRecord const * , std::pair<double,double> >;
static
reimputed_tags_in_this_window_t
reimpute_tags_one_by_one   (   mvn:: SquareMatrix const & C
                                , mvn:: SquareMatrix const & C_inv
                                , mvn:: SquareMatrix const & C1e8_inv
                                , std:: vector<double> const & tag_zs_
                                , std:: vector<RefRecord const *> const & rrs
                                , int const N_ref
                                , int const number_of_effective_tests_in_C_nolambda
                                );

enum class which_direction_t { DIRECTION_SHOULD_BE_REVERSED
                             , NO_ALLELE_MATCH
                             , DIRECTION_AS_IS };
string to_string(which_direction_t dir) {
    switch(dir) {
        break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED: return "DIRECTION_SHOULD_BE_REVERSED";
        break; case which_direction_t:: NO_ALLELE_MATCH:              return "NO_ALLELE_MATCH";
        break; case which_direction_t:: DIRECTION_AS_IS:              return "DIRECTION_AS_IS";
    }
    assert(1==2);
    return ""; // should never get here
}
static
    which_direction_t decide_on_a_direction
    ( RefRecord                       const & //r
    , SNPiterator     const & //g
    );


} // namespace ssimp

static
void    set_appropriate_locale(ostream & stream) {
    // By default, apply the user's locale, for example the thousands separator ...
    stream.imbue(std::locale(""));

    /*  However, within our test scripts, we force apostrophe ['] as the thousands
     *  separator so that we can compare results consistently across different developers
     *
     *  More info: http://en.cppreference.com/w/cpp/locale/numpunct/thousands_sep
     */

    if(std:: getenv("FORCE_THOUSANDS_SEPARATOR") != nullptr) {
        auto force_thousands_separator = std::string(std:: getenv("FORCE_THOUSANDS_SEPARATOR"));
        if(force_thousands_separator.size() ==1) {
            struct custom_thousands_separator : std::numpunct<char> {
                char    m_custom_separator;

                custom_thousands_separator(char c) : m_custom_separator(c) {}

                char do_thousands_sep()   const { return m_custom_separator; }
                std::string do_grouping() const { return "\3"; } // groups of 3 digits
            };
            stream.imbue(   std::locale (   stream.getloc()
                                        ,   new custom_thousands_separator(force_thousands_separator.at(0))
                                        // Nobody likes 'new' nowadays, but apparently this won't
                                        // leak. From http://en.cppreference.com/w/cpp/locale/locale/locale:
                                        //
                                        //  "Overload 7 is typically called with its second argument, f, obtained directly from a
                                        //   new-expression: the locale is responsible for calling the matching delete from its
                                        //   own destructor."
                                        )
                        );
        }
    }
}

bool exitWithUsage() {
    std:: cerr << usage_text;
    exit(1);
    return false;
}

int main(int argc, char **argv) {
    if(argc==1) {
        exitWithUsage();
    }

    options:: read_in_all_command_line_options(argc, argv);

    // There should be exactly zero, or three, non-option arguments:
    //     ssimp:  --gwas GWASFILE --ref REFPANEL --out IMPUTATIONOUTPUT
    // or
    //     ssimp:  GWASFILE REFPANEL IMPUTATIONOUTPUT
    // In the next few lines, we check this and ensure the three values
    // are stored in the appropriate variables.

    if( options:: opt_non_options.size() == 0 ) {
        // the three should already have been initialized via
        // --gwas, --ref, and --out:
        (!options:: opt_gwas_filename.empty()) || exitWithUsage();
        (!options:: opt_raw_ref      .empty()) || exitWithUsage();
        (!options:: opt_out          .empty()) || exitWithUsage();
    }
    else {
        options:: opt_non_options.size() == 3 || exitWithUsage();
        ( options:: opt_gwas_filename.empty()) || exitWithUsage();
        ( options:: opt_raw_ref      .empty()) || exitWithUsage();
        ( options:: opt_out          .empty()) || exitWithUsage();
        options:: opt_gwas_filename = options:: opt_non_options.at(0);
        options:: opt_raw_ref       = options:: opt_non_options.at(1);
        options:: opt_out           = options:: opt_non_options.at(2);
    }

    // Next we deal with the fact that the ref panel argument is
    // very special. If it begins with '1KG/', then we take the following
    // steps:
    //     1) Either the string is *exactly* '1KG' and the user must
    //        specify a --sample.names filename containing the list of
    //        sample IDS
    //     2) or, the user must *not* specify --sample.names and instead
    //        the sample IDs are inferred from the refpanel string:
    //            1KG/EUR
    //            1KG/TSI,EAS
    //            1KG/male
    //        The list of strings after the '/', separated by commas,
    //        is looked up in the 1kg panel file. If any of those
    //        strings appears in any of the columns, then that
    //        sample is used.
    if(options:: opt_raw_ref == "1KG") {
        !options:: opt_sample_names.empty() || DIE("If refpanel is exactly '1KG', then you must specify a --sample.names filename specifying the individuals to use");
    }
    if(options:: opt_raw_ref.substr(0, 4) == "1KG/") {
         options:: opt_sample_names.empty() || DIE("If refpanel begins with '1KG/', e.g. '1KG/EUR', then you must *not* specify --sample.names");
         options:: opt_raw_ref == "1KG/EUR" || DIE("Currently, only '1KG/EUR' is supported in this context");
         options:: opt_raw_ref      = "/home/amcdaid/reference_panels/1KG/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
         options:: opt_sample_names = "/home/amcdaid/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR";
    }

    if(!options:: opt_sample_names.empty()) {
        options:: adjust_sample_names_if_it_is_magical(); // "panelfile.txt/fieldname/filterfield=filtervalue"
    }


    if(!options:: opt_log.empty()) {
        logging:: setup_the_console_logging();
    }

    options:: list_of_tasks_to_run_at_exit.push_back(
        [](){
            auto PROF_CHANGE_DIR_AT_THE_LAST_MINUTE = getenv("PROF_CHANGE_DIR_AT_THE_LAST_MINUTE"); // to control where 'gmon.out' goes
            if(PROF_CHANGE_DIR_AT_THE_LAST_MINUTE) {
                int ret = chdir(PROF_CHANGE_DIR_AT_THE_LAST_MINUTE);
                (void)ret; // we don't care about this. But we need it anyway regarding a `warn_unused_result` attribute
            }
        });
    // This call to 'atexit' must come after 'setup_the_console_logging', in order to avoid
    // trying to delete the temporary file twice.
    atexit( [](void){
            for (   auto lifo = options:: list_of_tasks_to_run_at_exit.rbegin()
                ;        lifo != options:: list_of_tasks_to_run_at_exit.rend()
                ;      ++lifo) {
                (*lifo)();
            }
        });

    set_appropriate_locale(cout);

    cout << std::setprecision(20);

#if 0
    // hack for hpc1. If --ref isn't specified, default it to 1kg
    if( options:: opt_raw_ref.empty() ) {
        options:: opt_raw_ref="/data/sgg/aaron/shared/ref_panels/1kg/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
        cout << '\n';
        cout
            << "Note: as you didn't specify --ref, it is defaulted to 1KG '" << options:: opt_raw_ref << "'.\n"
            ;
        // Also, in this case, default the sample names to 'super_pop=EUR'
        if( options:: opt_sample_names.empty()) {
            options:: opt_sample_names = "/data/sgg/aaron/shared/ref_panels/1kg/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR";
            options:: adjust_sample_names_if_it_is_magical();
            cout
            << "      Also, as you didn't specify --sample.names, it has been defaulted to [" << options:: opt_sample_names << "].\n"
            << "      You can adjust that for other populations, e.g. \n"
            << "          --sample.names /data/sgg/aaron/shared/ref_panels/1kg/integrated_call_samples_v3.20130502.ALL.panel/sample/pop=TSI\n"
            ;
        }
        cout << '\n';
    }
#endif

    if( options:: opt_raw_ref.empty() ||  options:: opt_gwas_filename.empty()) {
        exitWithUsage();
        //DIE("Should pass args.\n    Usage:   " + string(argv[0]) + " --ref REFERENCEVCF --gwas GWAS --lambda 0.0 --window.width 1000000 --flanking.width 250000");
    }

    if( !options:: opt_impute_snps.empty() ) {
        gz:: igzstream f(options:: opt_impute_snps.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << options:: opt_impute_snps << ']');
        string one_SNPname_or_chrpos;
        options:: opt_impute_snps_as_a_uset  && DIE("multiple --impute.snps files not yet supported");
        options:: opt_impute_snps_as_a_uset = std::make_unique< std::unordered_set<std::string> >();
        while(getline( f, one_SNPname_or_chrpos)) {
            options:: opt_impute_snps_as_a_uset->insert(one_SNPname_or_chrpos);
        }
        options:: opt_impute_snps_as_a_uset->empty() && DIE("Nothing in --impute.snps file? [" << options:: opt_impute_snps << "]");
    }
    if( !options:: opt_tag_snps.empty() ) {
        gz:: igzstream f(options:: opt_tag_snps.c_str());
        (f.rdbuf() && f.rdbuf()->is_open()) || DIE("Can't find file [" << options:: opt_tag_snps << ']');
        string one_SNPname_or_chrpos;
        options:: opt_tag_snps_as_a_uset  && DIE("multiple --tag.snps files not yet supported");
        options:: opt_tag_snps_as_a_uset = std::make_unique< std::unordered_set<std::string> >();
        while(getline( f, one_SNPname_or_chrpos)) {
            options:: opt_tag_snps_as_a_uset->insert(one_SNPname_or_chrpos);
        }
        options:: opt_tag_snps_as_a_uset->empty() && DIE("Nothing in --tag.snps file? [" << options:: opt_tag_snps << "]");
    }

    if( options:: opt_out .empty()) { // --out wasn't specified - default to ${gwas}.ssimp.txt
        options:: opt_out = AMD_FORMATTED_STRING("{0}.ssimp.txt", options:: opt_gwas_filename);
    }
    if( options:: opt_log .empty()) { // --log wasn't specified - default to ${gwas}.log
        options:: opt_log = AMD_FORMATTED_STRING("{0}.log", options:: opt_gwas_filename);
    }

    // Finished adjusting the command line options

    // Now, do the imputation!

    if(!options:: opt_raw_ref.empty() && !options:: opt_gwas_filename.empty()) {
        auto gwas         = file_reading:: read_in_a_gwas_file(options:: opt_gwas_filename);

        // Go through regions, printing how many
        // SNPs there are in each region
        ssimp:: impute_all_the_regions(options:: opt_raw_ref, gwas);
    }
}

namespace ssimp{

static
chrpos  parse_chrpos_text_to_object     (string as_text, bool to_end_of_chromosome) {
    if  (   as_text.substr(0, 3) == "chr"
         || as_text.substr(0, 3) == "Chr"
         || as_text.substr(0, 3) == "CHR"
        )
        as_text = as_text.substr(3, string::npos);
    auto separated_by_colon = utils:: tokenize(as_text,':');
    chrpos position;
    switch(separated_by_colon.size()) {
        break; case 1: // no colon, just a chromosome
                position.chr = utils:: lexical_cast<int>(separated_by_colon.at(0));
                position.pos = to_end_of_chromosome
                                ?  std:: numeric_limits<int>::max()
                                :  std:: numeric_limits<int>::lowest() ;
        break; case 2:
                position.chr = utils:: lexical_cast<int>(separated_by_colon.at(0));
                position.pos = utils:: lexical_cast<int>(separated_by_colon.at(1));
        break; default:
                DIE("too many colons in [" << as_text << "]");
    }
    return position;
}
static
std:: pair<chrpos,chrpos> parse_range_as_pair_of_chrpos(string const &as_text) {
    auto separated_by_hyphen = utils:: tokenize(as_text, '-');
    switch(separated_by_hyphen.size()) {
        break; case 1: { // no hyphen. Just one position, or one entire chromosom
                auto lower = parse_chrpos_text_to_object(separated_by_hyphen.at(0), false);
                auto upper = parse_chrpos_text_to_object(separated_by_hyphen.at(0), true);
                return {lower,upper};
        }
        break; case 2: {
                auto lower = parse_chrpos_text_to_object(separated_by_hyphen.at(0), false);
                auto upper = parse_chrpos_text_to_object(separated_by_hyphen.at(1), true);
                return {lower,upper};
        }
        break; default:
                DIE("too many hyphens in [" << as_text << "]");
    }
    assert(0); // won't get here
    return { {-100,-100}, {100,1000} };
}

enum class enum_tag_or_impute_t { TAG, IMPUTE };
static
bool    test_if_skip(enum_tag_or_impute_t toi, RefRecord const &rr, int chrm) {

    // First, deal with --impute.snps or --tag.snps
    auto & up_snps = toi == enum_tag_or_impute_t::TAG ? options:: opt_tag_snps_as_a_uset : options:: opt_impute_snps_as_a_uset;
    if  (       up_snps ) {
        if (    up_snps->count( rr.ID ) == 0
             && up_snps->count( AMD_FORMATTED_STRING("chr{0}:{1}", chrm, rr.pos)) == 0
             ) {
            return true;
        }
    }

    // Next, deal with --impute.maf or --tag.maf
    auto & up_maf   = toi == enum_tag_or_impute_t::TAG
                    ? options:: opt_tag_maf
                    : options:: opt_impute_maf;
    if  ( rr.maf < up_maf ) {
            return true;
    }

    // --{impute,tag}.range
    auto & up_range = toi == enum_tag_or_impute_t::TAG
                    ? options:: opt_tag_range
                    : options:: opt_impute_range;
    if(!up_range.empty()) {
        auto up_range_parsed =  parse_range_as_pair_of_chrpos(up_range);
        using utils:: operator<<;
        if( chrpos{chrm, rr.pos} < up_range_parsed.first )
            return true;
        if( chrpos{chrm, rr.pos} > up_range_parsed.second )
            return true;
    }

    return false;
}
static
bool    test_if_skip(enum_tag_or_impute_t toi, chrpos l, chrpos u) {
    // 'l' and 'u' represent a window, or an entire chromosome.
    // We check here if --{impute,tag}.range is willing to skip
    // everything between l and u

    auto & up_range = toi == enum_tag_or_impute_t::TAG
                    ? options:: opt_tag_range
                    : options:: opt_impute_range;
    if(!up_range.empty()) {
        auto up_range_parsed =  parse_range_as_pair_of_chrpos(up_range);
        using utils:: operator<<;
        if( u < up_range_parsed.first )
            return true;
        if( l > up_range_parsed.second )
            return true;
    }

    return false;
}

static
void impute_all_the_regions(   string                                   filename_of_vcf
                             , file_reading:: GwasFileHandle_NONCONST   gwas
                             ) {
    unique_ptr<ofstream>   out_stream_for_imputations;
    ostream              * out_stream_ptr = nullptr;
    if(!options:: opt_out.empty()) {
        out_stream_for_imputations = make_unique<ofstream>( options:: opt_out. c_str() );
        (*out_stream_for_imputations) || DIE("Couldn't open the --out file [" + options:: opt_out + "]");
        out_stream_ptr = &*out_stream_for_imputations;
    }
    else {
        // If --out isn't specified, just print to stdout
        out_stream_ptr = & cout;
    }
    assert(out_stream_ptr);
    {   // print the header line
                    (*out_stream_ptr)
                        << "chr"
                        << '\t' << "pos"
                        << '\t' << "z_imp"
                        << '\t' << "source"
                        << '\t' << "SNP"
                        << '\t' << gwas->get_column_name_allele_ref() // copy column name from the GWAS input
                        << '\t' << gwas->get_column_name_allele_alt() // copy column name from the GWAS input
                        << '\t' << "maf"
                        << '\t' << "r2.pred"
                        << '\t' << "lambda"
                        << '\t' << "Z_reimputed"
                        << '\t' << "r2_reimputed"
                        << endl;
    }

    struct one_tag_data {   /* This is just to store data for --opt_tags_used_output for printing at the end,
                             * this is *not* used for actual imputation. */
        string      m_ID; // might be blank
        chrpos      m_chrpos;
        string      m_all_ref; // a1, other_allele , ...
        string      m_all_alt; // a2, effect_allele , ...
        double      m_z;
        double      m_N;

        bool operator< (one_tag_data const & other) const {
            if(m_chrpos     != other.m_chrpos   ) return m_chrpos   < other.m_chrpos    ;
            if(m_all_ref    != other.m_all_ref  ) return m_all_ref  < other.m_all_ref   ;
            if(m_all_alt    != other.m_all_alt  ) return m_all_alt  < other.m_all_alt   ;
            // The above three should be sufficient I think
            if(m_ID         != other.m_ID       ) return m_ID       < other.m_ID        ;
            return false;
        }
        bool operator==(one_tag_data const & other) const {
            return(
              (m_chrpos     == other.m_chrpos   )
            &&(m_all_ref    == other.m_all_ref  )
            &&(m_all_alt    == other.m_all_alt  )
            &&(m_ID         == other.m_ID       )
            );
        }
    };
    std:: vector<one_tag_data> tag_data_used;

    PP(options:: opt_window_width
      ,options:: opt_flanking_width
            );

    double N_max = -1;
    for(int i = 0; i<gwas->number_of_snps(); ++i ) {
        N_max = std::max(N_max, gwas->get_N(i));
    }

    int N_reference = -1; // to be updated (and printed) when we read in the first row of reference data
    int number_of_windows_seen_so_far_with_at_least_two_tags = 0; // useful to help decide when to do reimputation
    for(int chrm =  1; chrm <= 22; ++chrm) {
        cout.flush(); std::cerr.flush(); // helps with --log

        { // *optional* checks to help efficiency. Not necessary for correctness

        /* Skip entire chromosome? if --tags.range *or* --impute.range want
         * to skip this entire chromosome, then do it.
         */
            auto lowest_in_this_chromosome = chrpos{ chrm, std::numeric_limits<int>::lowest() };
            auto    max_in_this_chromosome = chrpos{ chrm, std::numeric_limits<int>::max()    };
            if  (   test_if_skip( enum_tag_or_impute_t:: IMPUTE
                        , lowest_in_this_chromosome, max_in_this_chromosome)
                 || test_if_skip( enum_tag_or_impute_t:: TAG
                        , lowest_in_this_chromosome, max_in_this_chromosome)
                )
                continue;
        }

        tbi:: read_vcf_with_tbi ref_vcf { filename_of_vcf, chrm };

        for(int w = 0; ; ++w ) {
            cout.flush(); std::cerr.flush(); // helps with --log

            int current_window_start = w     * options:: opt_window_width;
            int current_window_end   = (w+1) * options:: opt_window_width;

            /* Crucially, first we check if we have already gone past
             * the end of the current chromosome
             */
            {
                ref_vcf.set_region  (   chrpos{chrm,current_window_start - options:: opt_flanking_width}
                                    ,   chrpos{chrm,std::numeric_limits<int>::max()                    }
                                    );
                RefRecord rr;
                if(!ref_vcf.read_record_into_a_RefRecord(rr) ){
                    break; /* There is nothing in this window, *nor* later windows, hence
                            * we can break out of this chromosome
                            */
                }

            }

            // Next, check if we can skip the entire window due to --impute.range or --tag.range
            if  (   test_if_skip( enum_tag_or_impute_t:: IMPUTE
                        , chrpos{chrm,current_window_start - options:: opt_flanking_width}
                        , chrpos{chrm,current_window_end   + options:: opt_flanking_width} )
                 || test_if_skip( enum_tag_or_impute_t:: TAG
                        , chrpos{chrm,current_window_start}
                        , chrpos{chrm,current_window_end  } )
                )
                continue;

            /*
             * For a given window, there are three "ranges" to consider:
             *  -   The range of GWAS SNPs in the "broad" window (i.e. including the flanking region)
             *  -   The range of reference SNPs in the "broad" window.
             *  -   The range of reference SNPs in the "narrow" window (i.e. without the flanking region)
             *
             *  The intersection of the first two of these three ranges is the set of SNPs that
             *  are useful as tags. The third range is the set of targets.
             *
             *  (For now, we're including every tag as a target too.
             *
             *  In the next few lines, we find the endpoints (begin and end) of each of these
             *  three ranges.
             */

            // First, load up all the reference data in the broad window.

            vector<RefRecord>   all_nearby_ref_data;
            {   // Read in all the reference panel data in this broad window, only bi-allelic SNPs.
                ref_vcf.set_region  (   chrpos{chrm,current_window_start - options:: opt_flanking_width} // inclusive
                                    ,   chrpos{chrm,current_window_end   + options:: opt_flanking_width} // exclusive
                                    );
                RefRecord rr;
                while(ref_vcf.read_record_into_a_RefRecord(rr)) {
                    if(N_reference == -1) {
                        N_reference = rr.z12.size();
                        PP(N_reference);
                    }
                    else {
                        N_reference == ssize(rr.z12) || DIE("Varying size of reference panel? "
                                << N_reference
                                << " vs "
                                << rr.z12.size()
                                );
                    }

                    if  (   test_if_skip( enum_tag_or_impute_t:: IMPUTE, rr , chrm)
                         && test_if_skip( enum_tag_or_impute_t:: TAG   , rr , chrm)
                        )
                    {
                        /* skip this SNP as it's neither useful as a tag nor a target */
                        continue;
                    }
                    all_nearby_ref_data.push_back(rr);

                }

                for(int o=0; o+1 < ssize(all_nearby_ref_data); ++o) { // check position is monotonic
                    assert(all_nearby_ref_data.at(o).pos <= all_nearby_ref_data.at(o+1).pos);
                }
            }

            if(all_nearby_ref_data.empty()) {
                // no tags
                continue; // to the next window in this chromosome
            }

            {   // Update the gwas data with position information from the ref panel
                // I should put this chunk of code into another function
                std:: unordered_map<string, int> map_SNPname_to_pos; // just for this broad window
                for(auto & rr : all_nearby_ref_data) {
                    switch(map_SNPname_to_pos.count(rr.ID)) {
                        break; case 0:
                            map_SNPname_to_pos[rr.ID] = rr.pos;
                        break; case 1:
                        {
                            if(map_SNPname_to_pos[rr.ID] != rr.pos) {
                                PP(rr.ID); // TODO: remove this. I think it will be ID='.'
                                map_SNPname_to_pos[rr.ID] = -1;
                            }
                        }
                    }
                }
                auto it = map_SNPname_to_pos.begin();
                while(it != map_SNPname_to_pos.end()) {
                    if(it->second == -1) {
                        it = map_SNPname_to_pos.erase(it);
                    }
                    else {
                        ++it;
                    }
                }

                for(auto & msp : map_SNPname_to_pos) {
                    assert(msp.second != -1); // TODO: remove these instead
                }

                // Now that we have the position data from ref panel, see if we can copy it into the GWAS data
                auto gwas_it = begin_from_file(gwas);
                auto const e_gwas =   end_from_file(gwas);
                for(; gwas_it != e_gwas; ++gwas_it) {
                    auto gwas_SNPname = gwas_it.get_SNPname();
                    if(                     map_SNPname_to_pos.count( gwas_SNPname )) {
                        auto pos_in_ref =   map_SNPname_to_pos      [ gwas_SNPname ];
                        if  (   gwas_it.get_chrpos().chr == chrm
                             && gwas_it.get_chrpos().pos == pos_in_ref) {
                            // this fine, nothing to change, the two data sources agree
                        }
                        else if (   gwas_it.get_chrpos().chr == -1
                                 && gwas_it.get_chrpos().pos == -1) {
                            // just copy it in
                                gwas_it.set_chrpos(chrpos{chrm,pos_in_ref});
                        }
                        else
                            DIE("disagreement on position");
                    }
                }
            }
            //PP(__LINE__, utils:: ELAPSED());
            // Now, thanks to the block just completed, we have all the positions in the GWAS
            // that are relevant for this window.
            gwas->sort_my_entries(); // Sort them by position
            {   /* What follows is another optional speed optimization:
                 * If we now know all the positions of the tags, we may be able
                 * to skip windows (on the current chromome)
                 */
                if(gwas->get_chrpos(0) != chrpos{-1,-1}) { // every GWAS SNP has known position now
                    // check if the last tag on this chromosome is already too far back
                    auto x = std:: lower_bound( begin_from_file(gwas), end_from_file(gwas), chrpos{ chrm, std::numeric_limits<int>::lowest() });
                    auto last_tag_on_this_chromosome = chrpos{ chrm, std::numeric_limits<int>::lowest() };
                    while(x != end_from_file(gwas) && x.get_chrpos().chr == chrm) {
                        auto cur = x.get_chrpos();
                        assert( cur >= last_tag_on_this_chromosome );
                        last_tag_on_this_chromosome = cur;
                        ++x;
                    }
                    assert(chrm == last_tag_on_this_chromosome.chr);
                    if(last_tag_on_this_chromosome.pos < current_window_start - options:: opt_flanking_width) {
                        //PPe("break out", last_tag_on_this_chromosome.pos , current_window_start - options:: opt_flanking_width);
                        break;
                    }
                }
            }

            // Of the next four lines, only the two are interesting (to find the range of GWAS snps for this window)
            auto const b_gwas = begin_from_file(gwas);
            auto const e_gwas =   end_from_file(gwas);
            auto w_gwas_begin = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_start - options:: opt_flanking_width});
            auto w_gwas_end   = std:: lower_bound(b_gwas, e_gwas, chrpos{chrm,current_window_end   + options:: opt_flanking_width});

            // Next, find tags
            vector<double>                          tag_zs_;
            vector<RefRecord const *              > tag_its_;
            vector<double>                          tag_Ns;
            {   // Look for tags in the broad window.
                // This means moving monotonically through 'all_nearby_ref_data' and
                // the gwas data in parallel, and identifying suitable 'pairs'
                for(auto tag_candidate = range:: range_from_begin_end( w_gwas_begin, w_gwas_end )
                        ; ! tag_candidate.empty()
                        ;   tag_candidate.advance()
                        ) {
                    auto crps = tag_candidate.current_it().get_chrpos();
                    // Find the ref panel entries at the exact same chr:pos
                    auto ref_candidates = range:: range_from_begin_end(
                             std:: lower_bound( all_nearby_ref_data.begin() , all_nearby_ref_data.end(), chrpos(crps) )
                            ,std:: upper_bound( all_nearby_ref_data.begin() , all_nearby_ref_data.end(), chrpos(crps) )
                            );
                    vector<double> one_tag_zs;
                    vector<RefRecord const *> one_tag_its;
                    vector<int              > one_tag_Ns;
                    for(; ! ref_candidates.empty(); ref_candidates.advance()) {
                        auto & current_ref = ref_candidates.front_ref();
                        assert( current_ref.pos == crps.pos );

                        // apply --tag.{snps,maf,range}
                        if( test_if_skip( enum_tag_or_impute_t:: TAG, current_ref , chrm) ) {
                            continue;
                        }

                        auto dir = decide_on_a_direction( current_ref
                                             , tag_candidate .current_it() );
                        switch(dir) {
                            break; case which_direction_t:: DIRECTION_SHOULD_BE_REVERSED:
                                one_tag_zs.push_back( -tag_candidate.current_it().get_z() );
                                one_tag_its.push_back( &current_ref        );
                                one_tag_Ns.push_back(  tag_candidate.current_it().get_N() );
                            break; case which_direction_t:: DIRECTION_AS_IS             :
                                one_tag_zs.push_back(  tag_candidate.current_it().get_z() );
                                one_tag_its.push_back( &current_ref        );
                                one_tag_Ns.push_back(  tag_candidate.current_it().get_N() );
                            break; case which_direction_t:: NO_ALLELE_MATCH             : ;
                        }
                    }
                    // Finally, we make the decision on whether to include this tag or not.
                    // A tag is useful if there is exactly one reference SNP that corresponds
                    // to it (where 'correspond' means that the positions and alleles match).
                    assert(one_tag_zs.size() == one_tag_its.size());
                    if(1==one_tag_zs.size()) {
                        for(auto z : one_tag_zs)
                            tag_zs_.push_back(z);
                        for(auto it : one_tag_its)
                            tag_its_.push_back(it);
                        for(auto it : one_tag_Ns)
                            tag_Ns.push_back(it);
                    }
                }
            }

            assert(tag_zs_.size() == tag_its_.size());

            int const number_of_tags = tag_zs_.size();
            if(number_of_tags == 0)
                continue; // to the next window

            // Now, find suitable targets - i.e. anything in the reference panel in the narrow window
            // But some SNPs will have to be dropped, depending on --impute.range and --impute.snps
            vector<RefRecord const*>                unk2_its;
            {
                auto w_ref_narrow_begin = std:: lower_bound (   all_nearby_ref_data.begin() ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_start});
                auto w_ref_narrow_end   = std:: lower_bound (   all_nearby_ref_data.begin() ,   all_nearby_ref_data.end()
                                                            ,   chrpos{chrm,current_window_end});
                for(auto it = w_ref_narrow_begin; it<w_ref_narrow_end; ++it) {

                    // apply --impute.{snps,maf,range}
                    if( test_if_skip( enum_tag_or_impute_t:: IMPUTE, *it , chrm) ) {
                        continue;
                    }

                    auto const & z12_for_this_SNP = it->z12;
                    auto z12_minmax = minmax_element(z12_for_this_SNP.begin(), z12_for_this_SNP.end());
                    if  (*z12_minmax.first == *z12_minmax.second){
                        continue;   // no variation in the allele counts for this SNP within
                                    //the ref panel, therefore useless for imputation
                                    //(Actually, I think this might never happen, as the DISCARD rule will
                                    // take this into account. Anyway, no harm to leave this in.
                    }
                    unk2_its    .push_back( &*it       );
                }
            }
            //PP(__LINE__, utils:: ELAPSED());


            // Next few lines are kind of boring, just printing some statistics.

            // We have at least one SNP here, so let's print some numbers about this region
            auto number_of_snps_in_the_gwas_in_this_region      = w_gwas_end - w_gwas_begin;
            int  number_of_all_targets                          = unk2_its.size();

            if(number_of_all_targets == 0)
                continue;


            // We now have at least one tag, and at least one target, so we can proceed

            if  (number_of_tags > 1)
                ++number_of_windows_seen_so_far_with_at_least_two_tags;

            // TODO: I still am including all the tags among the targets - change this?

            cout
                << '\n'
                << "chrm" << chrm
                << "\t   " << current_window_start << '-' << current_window_end
                << '\n';

            cout << setw(8) << number_of_snps_in_the_gwas_in_this_region      << " # GWAS     SNPs in this window (with "<<options:: opt_flanking_width<<" flanking)\n";
            cout << setw(8) << number_of_tags                                 << " # SNPs in both (i.e. useful as tags)\n";
            cout << setw(8) << number_of_all_targets                          << " # target SNPs (anything in narrow window, will include some tags)\n";

            // Next, actually look up all the relevant SNPs within the reference panel
            vector<vector<int> const *> genotypes_for_the_tags  = from:: vector(tag_its_) |view::map| [](RefRecord const *rrp) { return &rrp->z12; } |action::collect;
            vector<vector<int> const *> genotypes_for_the_unks  = from:: vector(unk2_its) |view::map| [](RefRecord const *rrp) { return &rrp->z12; } |action::collect;

            assert(number_of_tags        == utils:: ssize(genotypes_for_the_tags));
            assert(number_of_all_targets == utils:: ssize(genotypes_for_the_unks));

            static_assert( std:: is_same< vector<vector<int> const *> , decltype(genotypes_for_the_tags) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation
            static_assert( std:: is_same< vector<vector<int> const *> , decltype(genotypes_for_the_unks) >{} ,""); // ints, not doubles, hence gsl_stats_int_correlation

            assert(!genotypes_for_the_tags.empty());
            assert(!genotypes_for_the_unks.empty());

            int const N_ref = genotypes_for_the_tags.at(0)->size(); // the number of individuals
            assert(N_ref > 0);
            assert(N_ref == N_reference); // TODO: redundant to have two variables like this

            //PP(__LINE__, utils:: ELAPSED());

            /*
             * Next few lines do a lot. The compute correlation, applying lambda regularization, and do imputation:
             */
            double const lambda = [&]() {
                if(options:: opt_lambda == "2/sqrt(n)")
                    return 2.0 / sqrt(N_reference);
                return utils:: lexical_cast<double>(options:: opt_lambda);
            }();

            assert(number_of_all_targets == ssize(unk2_its));

            // Next, record the real z of the tags. Just so we can print it in the output instead of trying to impute it
            std:: unordered_map<RefRecord const *, double> map_of_ref_records_of_tags;
            range:: zip_val( from:: vector(tag_its_), from:: vector(tag_zs_))
            |action::unzip_foreach|
            [&](auto && tag_refrecord, auto && z) {
                (void)z;
                map_of_ref_records_of_tags[tag_refrecord] = z;
            };
            assert(map_of_ref_records_of_tags.size() == tag_its_.size());

            // Compute the correlation matrices
            mvn:: SquareMatrix  C_nolambda  = make_C_tag_tag_matrix(genotypes_for_the_tags);
            mvn:: Matrix        c_nolambda  = make_c_unkn_tags_matrix( genotypes_for_the_tags
                                                         , genotypes_for_the_unks
                                                         , tag_its_
                                                         , unk2_its
                                                         );

            // Apply lambda
            auto                C_lambda        = apply_lambda_square(C_nolambda, lambda);
            auto                C_1e8lambda     = apply_lambda_square(C_nolambda, 1e-8);
            auto                c_lambda        = apply_lambda_rect(c_nolambda, unk2_its, tag_its_, lambda);
            auto                c_1e8lambda     = apply_lambda_rect(c_nolambda, unk2_its, tag_its_, 1e-8);


            if(options:: opt_missingness != options:: opt_missingness_t:: NAIVE) {
                for(double n : tag_Ns) {
                    n > 1 || DIE("--missingness requires n>1. [" << n << "]");
                    assert( n <= N_max );
                }
            }

            bool do_reimputation_in_this_window = false;
            if  (   number_of_tags > 1
                 && (   number_of_windows_seen_so_far_with_at_least_two_tags == 1 // this is the first such window
                     ||  options:: opt_reimpute_tags // if true, do reimputation in all such windows, not just the first one
                    )
                ) {
                do_reimputation_in_this_window = true;
            }

            // First, apply the missingness policy. Two matrices must be looped over
            std:: function<double(double,double,double)> current_missingness_policy;
            switch(options:: opt_missingness) {
                break;  case options:: opt_missingness_t:: NAIVE:
                    current_missingness_policy = [](double , double , double) -> double { return  1.0; }; // don't adjust C - just ignoring missingness for now
                break;  case options:: opt_missingness_t:: DEPENDENCY_MAXIMUM:
                    current_missingness_policy = [](double Nbig, double Nsml, double /*Nmax*/) -> double { assert(Nbig >= Nsml); return  std::sqrt(Nsml/Nbig); }; // assume maximum correlation in missingness
                break;  case options:: opt_missingness_t:: INDEPENDENCE:
                    current_missingness_policy = [](double Nbig, double Nsml, double Nmax) -> double { return  std::sqrt(Nsml)*std::sqrt(Nbig)/Nmax; };
            }
            for(int i=0; i<number_of_tags; ++i) {
                for(int j=0; j<number_of_tags; ++j) {
                    auto Ni = tag_Ns.at(i);
                    auto Nj = tag_Ns.at(j);
                    auto Nbig = std::max(Ni,Nj);
                    auto Nsml = std::min(Ni,Nj);
                    if(i!=j) {
                        auto delta_ij = current_missingness_policy(Nbig, Nsml, N_max); // *always* pass the big one in first
                        C_lambda.set   (i,j,  C_lambda   (i,j) * delta_ij);
                        C_1e8lambda.set(i,j,  C_1e8lambda(i,j) * delta_ij);
                    }
                }
            }
            assert(c_lambda.size1() == (size_t)number_of_all_targets); // number of rows
            assert(c_lambda.size2() == (size_t)number_of_tags); // number of columns
            for(int i=0; i<number_of_all_targets; ++i) {
                for(int j=0; j<number_of_tags; ++j) {
                    auto Nj = tag_Ns.at(j);
                    assert( Nj <= N_max );

                    auto delta_ij = current_missingness_policy(N_max, Nj, N_max); // *always* pass the big one in first
                    assert(delta_ij > 0.0);
                    assert(delta_ij <= 1.0);
                    c_lambda.set   (i,j,  c_lambda   (i,j) * delta_ij);
                    c_1e8lambda.set(i,j,  c_1e8lambda(i,j) * delta_ij);
                }
            }
            // missingness has now been applied.

            // Compute the imputations
            auto c_Cinv_zs = mvn:: multiply_matrix_by_colvec_giving_colvec(c_lambda, solve_a_matrix (C_lambda, mvn:: make_VecCol(tag_zs_)));

            int number_of_effective_tests_in_C_nolambda = compute_number_of_effective_tests_in_C_nolambda(C_1e8lambda);

            // if necessary, do reimputation
            // TODO: what is the r2 for the reimputation? I'm being 'naive' for now
            reimputed_tags_in_this_window_t reimputed_tags_in_this_window;
            if  (do_reimputation_in_this_window) {
                reimputed_tags_in_this_window = reimpute_tags_one_by_one(C_lambda, invert_a_matrix(C_lambda), invert_a_matrix(C_1e8lambda), tag_zs_, tag_its_, N_ref, number_of_effective_tests_in_C_nolambda);
            }
            assert(ssize(reimputed_tags_in_this_window) == 0
                || ssize(reimputed_tags_in_this_window) == number_of_tags);


            // Next few lines are for the imputation quality
            auto C1e8inv_c1e8   =     mvn:: multiply_NoTrans_Trans( invert_a_matrix(C_1e8lambda)  , c_1e8lambda);
            vector<double> imp_quals_corrected; // using the 'number_of_effective_tests' thing

            for(int i=0; i<number_of_all_targets; ++i) {
                mvn:: Matrix to_store_one_imputation_quality(1,1); // a one-by-one-matrix
                auto imp_qual = [&](){ // multiply one vector by another, to directly get
                                        // the diagonal of the imputation quality matrix.
                                        // This should be quicker than fully computing
                                        // c' * inv(C) * c
                    auto lhs = gsl_matrix_const_submatrix(   c_1e8lambda.get(), i, 0, 1, number_of_tags);
                    auto rhs = gsl_matrix_const_submatrix( C1e8inv_c1e8 .get(), 0, i, number_of_tags, 1);
                    // TODO: Maybe move the next line into mvn.{hh,cc}?
                    const int res_0 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &lhs.matrix, &rhs.matrix, 0, to_store_one_imputation_quality.get());
                    assert(res_0 == 0);
                    auto simple_imp_qual = to_store_one_imputation_quality(0,0);
                    //PPe(N_ref, simple_imp_qual);
                    return adjusted_imputation_quality(simple_imp_qual, N_ref, number_of_effective_tests_in_C_nolambda);
                }();
                imp_quals_corrected.push_back(imp_qual);
            }
            assert(number_of_all_targets == ssize(imp_quals_corrected));

            // Finally, print out everything to the --out file
            for(int i=0; i<number_of_all_targets; ++i) {
                    auto && target = unk2_its.at(i);
                    auto pos  = target->pos;
                    auto SNPname = target->ID;

                    auto imp_qual = imp_quals_corrected.at(i);

                    auto z_imp = c_Cinv_zs(i);
                    double Z_reimputed = std::nan("");
                    double r2_reimputed = std::nan("");

                    // Is this target actually a tag SNP?
                    bool const was_in_the_GWAS = map_of_ref_records_of_tags.count(target);

                    if(was_in_the_GWAS) { // don't print the imputation, just copy straight from the GWAS
                        auto GWAS_z = map_of_ref_records_of_tags.at(target);
                        // assert(std::abs(imp_qual-1.0) < 1e-5);
                        // assert(std::abs(GWAS_z-z_imp) < 1e-5); // TODO: breaks with 1KG/EUR + UKB-both - I should investigate this further
                        imp_qual = 1.0;
                        z_imp = GWAS_z;

                    }
                    if(was_in_the_GWAS && reimputed_tags_in_this_window.size() > 0) {
                        assert(reimputed_tags_in_this_window.count(target) == 1);
                        Z_reimputed  = reimputed_tags_in_this_window.at(target).first;
                        r2_reimputed = reimputed_tags_in_this_window.at(target).second;
                    }

                    (*out_stream_ptr)
                                << chrm
                        << '\t' << pos
                        << '\t' << z_imp
                        << '\t' << (was_in_the_GWAS ? "GWAS" : "SSIMP")
                        << '\t' << SNPname
                        << '\t' << target->ref
                        << '\t' << target->alt
                        << '\t' << target->maf
                        << '\t' << imp_qual
                        << '\t' << lambda
                        << '\t' << (std::isnan( Z_reimputed) ? "" : AMD_FORMATTED_STRING("{0}",  Z_reimputed))
                        << '\t' << (std::isnan(r2_reimputed) ? "" : AMD_FORMATTED_STRING("{0}", r2_reimputed))
                        << endl;
            }

            // Finally, store every tag used. We will print it all at the end to
            // another file (--tags.used.output)

            if(!options:: opt_tags_used_output.empty()) {
                assert(tag_its_.size() == tag_zs_.size());
                range:: zip_val ( from:: vector(tag_its_)
                                , from:: vector(tag_zs_)
                                , from:: vector(tag_Ns)
                                )
                |action::unzip_foreach|
                [&](auto && tag_refrecord, auto z, auto N) {
                    one_tag_data otd {
                        tag_refrecord->ID
                            , {chrm, tag_refrecord->pos}
                        , tag_refrecord->ref
                        , tag_refrecord->alt
                        , z
                        , N
                    };
                    tag_data_used.push_back(otd);
                };
            }
        } // loop over windows
    } // loop over chromosomes

    if(!options:: opt_tags_used_output.empty()) {
        ofstream opts_tags_used_output_file(options:: opt_tags_used_output); // to store the tags used, if requested by --opt_tags_used_output
        opts_tags_used_output_file || DIE("Couldn't create --tags.used.output [" << options:: opt_tags_used_output << "] for output");

        // Remove duplicates from the 'tag_data_used' structure, as a tag can be a tag
        // in multiple windows.
        sort(tag_data_used.begin(), tag_data_used.end());
        tag_data_used.erase(
                unique(tag_data_used.begin(), tag_data_used.end())
            ,   tag_data_used.end());

        opts_tags_used_output_file << std::setprecision(53);
        opts_tags_used_output_file
                    << "ID"
            << '\t' << "Chr"
            << '\t' << "Pos"
            << '\t' << gwas->get_column_name_allele_ref() // copy column name from the GWAS input
            << '\t' << gwas->get_column_name_allele_alt() // copy column name from the GWAS input
            << '\t' << "z"
            << '\t' << "N"
            << '\n'
            ;
        for(auto && otd : tag_data_used) {
            opts_tags_used_output_file
                    << otd.m_ID
            << '\t' << otd.m_chrpos.chr
            << '\t' << otd.m_chrpos.pos
            << '\t' << otd.m_all_ref
            << '\t' << otd.m_all_alt
            << '\t' << otd.m_z
            << '\t' << otd.m_N
            << '\n';
        }
    }
}

static
mvn:: SquareMatrix
make_C_tag_tag_matrix( vector<vector<int> const *>      const & genotypes_for_the_tags) {
    int const number_of_tags = genotypes_for_the_tags.size();
    int const N_ref = genotypes_for_the_tags.at(0)->size();
    assert(N_ref > 0);

    mvn:: SquareMatrix C (number_of_tags);
    for(int k=0; k<number_of_tags; ++k) {
        for(int l=0; l<number_of_tags; ++l) {
            assert(N_ref        == utils:: ssize(*genotypes_for_the_tags.at(k)));
            assert(N_ref        == utils:: ssize(*genotypes_for_the_tags.at(l)));
            double c_kl = gsl_stats_int_correlation( &genotypes_for_the_tags.at(k)->front(), 1
                                                   , &genotypes_for_the_tags.at(l)->front(), 1
                                                   , N_ref );
            if(c_kl > (1.0-1e-10)) { // sometimes it sneaks above one, don't really know how
                assert(c_kl-1.0 < 1e-5);
                assert(c_kl-1.0 >-1e-5);
                c_kl = 1.0;
            }
            if(c_kl < (-1.0+1e-10)) { // in case it goes below -1.0 also
                assert(c_kl+1.0 < 1e-5);
                assert(c_kl+1.0 >-1e-5);
                c_kl =-1.0;
            }
            assert(c_kl >= -1.0);
            assert(c_kl <=  1.0);
            if(k==l) {
                assert(c_kl == 1.0);
            }
            else {
            }
            C.set(k,l,c_kl);
        }
    }
    return C;
}
static
mvn:: Matrix make_c_unkn_tags_matrix
        ( vector<vector<int> const *>      const & genotypes_for_the_tags
        , vector<vector<int> const *>      const & genotypes_for_the_unks
        , vector<RefRecord const *              > const & tag_its
        , vector<RefRecord const *              > const & unk_its
        ) {
    int const number_of_tags = genotypes_for_the_tags.size();
    int const number_of_all_targets = genotypes_for_the_unks.size();
    assert(number_of_tags        == ssize(tag_its));
    assert(number_of_all_targets == ssize(unk_its));
    int const N_ref = genotypes_for_the_tags.at(0)->size();
    assert(N_ref > 0);

    mvn:: Matrix      c(number_of_all_targets, number_of_tags);

    for(auto tags : zip_val ( range:: ints(number_of_tags)
                            , range:: range_from_begin_end(tag_its)
                            , range:: range_from_begin_end(genotypes_for_the_tags)
                )) {
        int     k           = std::get<0>(tags);
        auto    tag_its_k   = std::get<1>(tags);
        auto &  calls_at_k  =*std::get<2>(tags);

        assert(&calls_at_k ==  genotypes_for_the_tags.at(k)); // double check that we're getting it by-reference, i.e via the .get() on the previous line

        assert(tag_its_k == tag_its.at(k));

        zip_val ( range:: ints(number_of_all_targets)
                , range:: range_from_begin_end(unk_its)
                , range:: range_from_begin_end(genotypes_for_the_unks) | view:: ref_wraps
                )
        |action:: unzip_foreach|
        [&] (   int
                    u
            ,   RefRecord const *
                    unk_its_u
            ,   vector<int> const *
                    calls_at_u_rw
            ) -> void {
                vector<int> const &  calls_at_u  = *calls_at_u_rw;
                assert(&calls_at_u ==  genotypes_for_the_unks.at(u));

                assert(N_ref        == utils:: ssize(calls_at_k));
                assert(N_ref        == utils:: ssize(calls_at_u));
                double c_ku = gsl_stats_int_correlation( &calls_at_k.at(0), 1
                                                       , &calls_at_u.at(0), 1
                                                       , N_ref );
                if(c_ku > (1.0-1e-10)) {
                    assert(c_ku-1.0 < 1e-5);
                    assert(c_ku-1.0 >-1e-5);
                    c_ku = 1.0;
                }
                if(c_ku < (-1.0+1e-10)) {
                    assert(c_ku+1.0 < 1e-5);
                    assert(c_ku+1.0 >-1e-5);
                    c_ku = -1.0;
                }

                if(tag_its_k == unk_its_u) { // identical pointer value, i.e. same reference entry
                    assert(c_ku ==  1.0);
                }
                else {
                    assert(c_ku >= -1.0);
                    assert(c_ku <=  1.0);
                }
                c.set(u,k,c_ku);
         };
    }
    return c;
};
static
mvn::SquareMatrix apply_lambda_square (mvn:: SquareMatrix copy, double lambda) {
    int s = copy.size();
    for(int i=0; i<s; ++i) {
        for(int j=0; j<s; ++j) {
            if(i!=j)
                copy.set(i,j,  copy(i,j)*(1.0-lambda));
        }
    }
    return copy;
}
static
mvn:: Matrix apply_lambda_rect   (   mvn:: Matrix copy
                        ,   vector<RefRecord const *> const & unk2_its
                        ,   vector<RefRecord const *> const & tag_its_
                        , double lambda) {
    int number_of_all_targets = unk2_its.size();
    int number_of_tags        = tag_its_.size();
    assert(copy.size1() == (size_t)number_of_all_targets);
    assert(copy.size2() == (size_t)number_of_tags);
    for(int i=0; i<number_of_all_targets; ++i) {
        for(int j=0; j<number_of_tags; ++j) {
            if(tag_its_.at(j) == unk2_its.at(i))
                assert(copy(i,j) == 1.0);
            else {
                copy.set(i,j, copy(i,j)*(1.0-lambda));
            }
        }
    }
    return copy;
}
static
int compute_number_of_effective_tests_in_C_nolambda(mvn:: SquareMatrix const & C_smalllambda) {
    int number_of_tags = C_smalllambda.size();
    vector<double> eigenvalues;

    mvn:: VecCol eig_vals (number_of_tags);
    // compute the eigenvalues of C. TODO: C_lambda or C_smalllambda?
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc(number_of_tags);
    auto A = C_smalllambda; // copy the matrix, so that it can be destroyed by gsl_eigen_symm
    int ret = gsl_eigen_symm (A.get(), eig_vals.get(), w);
    assert(ret==0);
    gsl_eigen_symm_free(w);
    for(int i = 0; i< utils::ssize(eig_vals); ++i)
        eigenvalues.push_back( eig_vals(i) );

    sort(eigenvalues.rbegin(), eigenvalues.rend());
    //using utils:: operator<<;
    //PPe(eigenvalues);
    auto sum_of_eigvals = std:: accumulate(eigenvalues.begin(), eigenvalues.end(), 0.0);
    assert( std::abs(sum_of_eigvals-number_of_tags) < 1e-5 );
    double running_total_of_eigenvalues = 0.0;

    for(int i = 0; i< utils::ssize(eigenvalues); ++i) {
        running_total_of_eigenvalues += eigenvalues.at(i);
        //PPe(running_total_of_eigenvalues);
        if( running_total_of_eigenvalues >= 0.995 * sum_of_eigvals ) {
            return i+1;
        }
    }
    assert(0);
    return -1; // never going to get here
}
static
double  adjusted_imputation_quality (double iq_1e8_simple, int N_ref, int number_of_effective_tests) {
    /* With a small change to the method which is in the paper currently.
     * As the number_of_effective_tests increases, the resulting (estimate of)
     * imputation quality decreases.
     *
     * However, a very large number_of_effective_tests could give extreme values, much larger than 1.0.
     * I've seen +400.0 and -400.0
     *
     * The IQ is zero when:
     *  1.0 = (1.0-iq_1e8_simple)*(N_ref - 1.0)/(N_ref - number_of_effective_tests - 1.0);
     *  (N_ref - number_of_effective_tests - 1.0) =             (1.0-iq_1e8_simple)*(N_ref - 1.0);
     *         - number_of_effective_tests        = 1 - N_ref + (1.0-iq_1e8_simple)*(N_ref - 1.0);
     *           number_of_effective_tests        = N_ref - 1 - (1.0-iq_1e8_simple)*(N_ref - 1.0);
     *           number_of_effective_tests        = N_ref - 1 + (iq_1e8_simple-1.0)*(N_ref - 1.0);
     *           number_of_effective_tests        =           (1+iq_1e8_simple-1.0)*(N_ref - 1.0);
     *           number_of_effective_tests        =              iq_1e8_simple     *(N_ref - 1.0);
     */
    if(number_of_effective_tests  >= iq_1e8_simple*(N_ref - 1.0)) {
        // A little over this threshold results in negative IQ. But if
        // it exceeds N_ref, then it becomes positive again!
        // In either case though, it should (Aaron thinks) return 0.0
        return 0;
    }
    return 1.0 - (1.0-iq_1e8_simple)*(N_ref - 1.0)/(N_ref - number_of_effective_tests - 1.0);
}

static
    which_direction_t decide_on_a_direction
    ( RefRecord                       const & r
    , SNPiterator     const & g
    ) {
        auto const rp_ref = r.ref;
        auto const rp_alt = r.alt;
        auto const gw_ref = g.get_allele_ref();
        auto const gw_alt = g.get_allele_alt();

        //PPe(r.pos ,rp_ref ,rp_alt ,gw_ref ,gw_alt);
        assert(rp_ref != rp_alt);
        assert(gw_ref != gw_alt);

        if(1)
        { // check for palindromes - we can't be sure of their direction
            if  (   rp_ref.find_first_of("AT") != string:: npos
                 && rp_alt.find_first_of("AT") != string:: npos) { // palindromic - can't really use it. TODO: Should this be a command line option?
                return which_direction_t:: NO_ALLELE_MATCH;
            }
            if  (   rp_ref.find_first_of("GC") != string:: npos
                 && rp_alt.find_first_of("GC") != string:: npos) { // palindromic - can't really use it. TODO: Should this be a command line option?
                return which_direction_t:: NO_ALLELE_MATCH;
            }
            if  (   gw_ref.find_first_of("AT") != string:: npos
                 && gw_alt.find_first_of("AT") != string:: npos) { // palindromic - can't really use it. TODO: Should this be a command line option?
                return which_direction_t:: NO_ALLELE_MATCH;
            }
            if  (   gw_ref.find_first_of("GC") != string:: npos
                 && gw_alt.find_first_of("GC") != string:: npos) { // palindromic - can't really use it. TODO: Should this be a command line option?
                return which_direction_t:: NO_ALLELE_MATCH;
            }
        }

        if( rp_ref == gw_ref
         && rp_alt == gw_alt)
            return which_direction_t:: DIRECTION_AS_IS;

        if( rp_ref == gw_alt
         && rp_alt == gw_ref)
            return which_direction_t:: DIRECTION_SHOULD_BE_REVERSED;

        return which_direction_t:: NO_ALLELE_MATCH;
    }

static
reimputed_tags_in_this_window_t
reimpute_tags_one_by_one   (   mvn:: SquareMatrix const & C
                                , mvn:: SquareMatrix const & C_inv
                                , mvn:: SquareMatrix const & C1e8_inv
                                , std:: vector<double> const & zs
                                , std:: vector<RefRecord const *> const & rrs
                                , int const N_ref
                                , int const number_of_effective_tests_in_C_nolambda
                                ) {
    assert(C.size() == zs.size());
    assert(C.size() == C_inv.size());
    assert(C.size() == rrs.size());
    int const M = C.size();
    assert(M>1);

    mvn:: VecCol vec(M);
    std:: vector<double> reimputed_tags_z;
    std:: vector<double> reimputed_tags_r2;
    reimputed_tags_in_this_window_t map_of_ref_records_of_reimputation_results;
    for(int m=0; m<M; ++m) {
        auto z_real = zs.at(m);
        auto * rrp = rrs.at(m);
        auto z_fast = [&]() {
            for(int n=0; n<M; ++n) {
                vec.set(n, C_inv(n,m));
            }
            vec.set(m,0); // to ensure the known z is ignored

            auto x = multiply_rowvec_by_colvec_giving_scalar( vec, mvn:: make_VecCol(zs) )(0) ;
            auto y =C_inv(m,m);

            auto reimp = -x/y;
            auto iq_lambda_simple = C(m,m)-1.0/y; // we won't actually report this any more

            // to do the '#effective tests' calculation of IQ, we need to start with the
            // simple IQ, *but* with negligible lamdba

            assert(C(m,m) == 1);
            auto y_1e8 =C1e8_inv(m,m);
            auto iq_1e8_simple = 1-1.0/y_1e8; // we won't actually report this any more

            auto iq_1e8_effective = adjusted_imputation_quality(iq_1e8_simple, N_ref, number_of_effective_tests_in_C_nolambda);

            (void)iq_lambda_simple;
            (void)iq_1e8_effective;

            return std:: make_pair(reimp
                    //, iq_lambda_simple
                    , iq_1e8_effective
                    );
        }();
        reimputed_tags_z .push_back(z_fast.first);
        reimputed_tags_r2.push_back(z_fast.second);
        map_of_ref_records_of_reimputation_results[ rrp ] = z_fast;
        (void)z_real;
#if 0
        PPe(M, m, z_real, z_fast.first, z_fast.second
                , rrp->maf
                , rrp->ID
                , rrp->pos
                , rrp->ref
                , rrp->alt
                );
        auto z_slow = [&]() {
            // I'll simply set the correlations to zero
            auto C_zeroed = C;
            for(int n=0; n<M; ++n) {
                if(n!=m) {
                    // the 'real' value has no correlation with the other tags
                    C_zeroed.set(n,m,0);
                    C_zeroed.set(m,n,0);
                }
                else
                    assert(C_zeroed(n,m)==1);
            }
            auto C_zeroed_inv = invert_a_matrix(C_zeroed);
            mvn:: VecCol vec(M);
            for(int n=0; n<M; ++n) {
                vec.set(n, C(n,m));
            }
            vec.set(m, 0); // and the target has no correlation with the 'real' value
            auto z_slow_imp
                =   multiply_rowvec_by_colvec_giving_scalar
                    (   vec
                    ,   multiply_matrix_by_colvec_giving_colvec
                        (   C_zeroed_inv
                        , mvn:: make_VecCol(zs)
                        )
                    );
            auto z_slow_IQ
                =   multiply_rowvec_by_colvec_giving_scalar
                    (   vec
                    ,   multiply_matrix_by_colvec_giving_colvec
                        (   C_zeroed_inv
                        , vec
                        )
                    );
            return std:: make_pair( z_slow_imp(0), z_slow_IQ(0) );
        }();
        assert(std:: fabs(z_slow.first - z_fast.first) < 1e-5);
        assert(std:: fabs(z_slow.second- z_fast.second)< 1e-5);
#endif
    }
    return map_of_ref_records_of_reimputation_results;
}

} // namespace ssimp
