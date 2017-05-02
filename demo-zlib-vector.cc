#include "zlib-vector.hh"

#include <string>
#include <cassert>

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

#include"../module-bits.and.pieces/DIE.hh"

using std:: vector;


/* Compress from file source to file dest until EOF on source.
   def() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_STREAM_ERROR if an invalid compression
   level is supplied, Z_VERSION_ERROR if the version of zlib.h and the
   version of the library linked do not match, or Z_ERRNO if there is
   an error reading or writing the files. */
static
vector<unsigned char> read_FILE_into_vector(FILE *source) {
    vector<unsigned char> c;

    constexpr int CHUNK = 16384;
    unsigned char tmp[CHUNK];

    do {
        size_t n = fread(tmp, 1, CHUNK, source);
        if (ferror(source)) {
            DIE("some problem reading from stdin");
        }
        if(n==0)
            break;
        c.insert( c.end(), tmp, tmp+n );
    } while(1);
    return c;
}

/* Decompress from file source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files. */

/* compress or decompress from stdin to stdout */
int main(int argc, char **argv)
{
    /* avoid end-of-line conversions */
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    /* do compression if no arguments */
    if (argc == 1) {
        auto input_data = read_FILE_into_vector(stdin);
        auto result = zlib_vector:: deflate(std::move(input_data));
        std:: cout.write    (   reinterpret_cast<char*>(result.data())
                            ,   result.size()
                            );
        return 0;
    }

    /* do decompression if -d specified */
    else if (argc == 2 && std::string(argv[1]) == "-d") {
        auto input_data = read_FILE_into_vector(stdin);
        auto result = zlib_vector:: inflate(std::move(input_data));
        std:: cout.write    (   reinterpret_cast<char*>(result.data())
                            ,   result.size()
                            );
        return 0;
    }

    /* otherwise, report usage */
    else {
        fputs("zpipe usage: zpipe [-d] < source > dest\n", stderr);
        return 1;
    }
}
