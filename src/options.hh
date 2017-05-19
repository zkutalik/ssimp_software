#include <string>

namespace options {

extern  std:: string            opt_raw_ref;
extern  std:: string            opt_gwas_filename;
extern  std:: string            opt_out;
extern  int                     opt_window_width;
extern  int                     opt_flanking_width;
extern  double                  opt_lambda;
extern  std:: string            opt_impute_range;
extern  std:: string            opt_impute_snps;

        void                    read_in_all_command_line_options(int argc, char **argv);

} // namespace options
