
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h> // for wait()
#include <sys/wait.h> // for wait()
#include <fcntl.h> // for F_GETFD

#include <fstream>

#include "bits.and.pieces/DIE.hh"
#include "bits.and.pieces/PP.hh"

#include "options.hh"

namespace logging {
void setup_the_console_logging() {
        std:: ofstream log(options:: opt_log);
        log || DIE("Couldn't open log file --log [" << options:: opt_log << "]");

        // Set up a child process, which will listen on two file descriptors,
        // forwarding everything to stdout and stderr while recording
        // everything in the --log file.
        //
        // The parent process (i.e. the main 'ssimp' process) but connect both
        // its outputs to this child process.

        pid_t mypipe_out[2]; // for stdout (cout)
        pid_t mypipe_err[2]; // for stderr (cerr)

        pipe(mypipe_out) && DIE("pipe() failed.");
        pipe(mypipe_err) && DIE("pipe() failed.");

        auto p = fork();
        p != -1 || DIE("Couldn't fork(). Too many processes");
        if(p== (pid_t)0) { // child
            close(mypipe_out[1]);
            close(mypipe_err[1]);

            bool keep_going;

            do {
                int nfds_for_select = std::max(mypipe_out[0], mypipe_err[0]) + 1;
                fd_set rfds;
                FD_ZERO(&rfds);
                FD_SET(mypipe_out[0], &rfds);
                FD_SET(mypipe_err[0], &rfds);

                select(nfds_for_select, &rfds, NULL, NULL, NULL);

                keep_going = false; // will be changed to 'true' if something is read on either fd

                char buf[1024];

                if(FD_ISSET(mypipe_out[0], &rfds)) {
                    int chars_read = read(mypipe_out[0], buf, sizeof(buf));
                    if(chars_read > 0) {
                        keep_going = true;
                        {int ret=write(1, buf, chars_read);(void)ret;}
                        log.write(buf, chars_read);
                        log.flush();
                    }
                }
                if(FD_ISSET(mypipe_err[0], &rfds)) {
                    int chars_read = read(mypipe_err[0], buf, sizeof(buf));
                    if(chars_read > 0) {
                        keep_going = true;
                        {int ret=write(2, buf, chars_read);(void)ret;}
                        log.write(buf, chars_read);
                        log.flush();
                    }
                }
            } while(keep_going);
            // At this point, both file descriptors have EOF

            log.close(); // just in case the following 'exit(0)' fails to flush and close
            exit(0); // exit the child
        }
        else { //parent
            close (mypipe_out[0]);
            close (mypipe_err[0]);

            // Connect my stdout to the child process
            dup2(mypipe_out[1], 1) == 1 || DIE("Couldn't dup2() into the pipe");
            // Connect my stderr to the child process
            dup2(mypipe_err[1], 2) == 2 || DIE("Couldn't dup2() into the pipe");

            // From here, the rest of the imputation program continues. All
            // further writes to 'cout' or to 'cerr' go through the child
            // process above, which logs everything to the file specified
            // by --log
        }
}
} //namespace logging
