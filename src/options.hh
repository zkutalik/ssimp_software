#include <string>
#include <unordered_set>
#include <functional>
#include <vector>
#include <set>
#include <memory> // for unique_ptr

namespace options {

extern  int                     opt_help /*= 0*/;
extern  std::vector<std:: string>           opt_non_options; // if there are two or three of these, copy them into opt_gwas_filename, opt_out, and opt_raw_ref

extern  std:: string            opt_raw_ref;
extern  std:: string            opt_gwas_filename;
extern  std:: string            opt_out;
extern  std:: string            opt_log; // copy of whatever is sent to the console.
extern  int                     opt_window_width;
extern  int                     opt_flanking_width;
extern  std:: string            opt_lambda;

extern  std:: string            opt_impute_range;
extern  std:: string            opt_impute_snps;
extern  std::unique_ptr<std::unordered_set<std::string>>    opt_impute_snps_as_a_uset;
extern  double                  opt_impute_maf/* =0.0*/; // target not imputed unless maf (in reference) is at least this.

        // The next few are like the --impute.* above, but applying to tags instead
extern  std:: string            opt_tag_range;
extern  std:: string            opt_tag_snps;
extern  std::unique_ptr<std::unordered_set<std::string>>    opt_tag_snps_as_a_uset;
extern  double                  opt_tag_maf/* =0.0*/;

extern  bool                    opt_reimpute_tags;
extern  std:: string            opt_tags_used_output;

extern  std:: string            opt_sample_names;
enum class opt_missingness_t { NAIVE, DEPENDENCY_MAXIMUM, INDEPENDENCE };
extern  opt_missingness_t       opt_missingness;


extern  std:: vector<std::function<void(void)>>    list_of_tasks_to_run_at_exit;

extern  int                     opt_download_build_db;
extern  int                     opt_download_1KG;
extern  std:: set<int>          opt_debug_build_chromosomes_to_load; // if empty, load all of them

        void                    read_in_all_command_line_options(int argc, char **argv);
        void                    adjust_sample_names_if_it_is_magical();


} // namespace options
