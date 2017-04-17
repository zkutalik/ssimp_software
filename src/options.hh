#include <string>

namespace options {

extern  std:: string            opt_raw_ref;
extern  std:: string            opt_gwas_filename;
extern  std:: string            opt_out;
extern  int                     opt_window_width;
extern  int                     opt_flanking_width;
extern  double                  opt_lambda;

        void                    read_in_all_command_line_options(int argc, char **argv);

} // namespace options
