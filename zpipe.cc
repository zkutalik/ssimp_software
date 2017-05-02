/* zpipe.c: example of proper use of zlib's inflate() and deflate()
   Not copyrighted -- provided to the public domain
   Version 1.4  11 December 2005  Mark Adler */

/* Version history:
   1.0  30 Oct 2004  First version
   1.1   8 Nov 2004  Add void casting for unused return values
                     Use switch statement for inflate() return values
   1.2   9 Nov 2004  Add assertions to document zlib guarantees
   1.3   6 Apr 2005  Remove incorrect assertion in inf()
   1.4  11 Dec 2005  Add hack to avoid MSDOS end-of-line conversions
                     Avoid some compiler warnings for input and output buffers
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "zlib.h"

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

#include<vector>
#include"../module-bits.and.pieces/PP.hh"
#include"../module-bits.and.pieces/DIE.hh"

using std:: vector;
using std:: cerr;

namespace exc {
    // Too lazy to make proper std::exception objects to throw. So
    // I'll just throw these instead:
    struct error_from_deflateInit_t {};
    struct error_from_inflateInit_t {
        decltype(Z_DATA_ERROR) m_error_code;
    };
}

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
    c.push_back('\0');
    cerr << "c.size()=" << c.size() << '\n';
    c.pop_back();
    return c;
}
vector<unsigned char> def(vector<unsigned char> src, int level)
{
    constexpr int CHUNKdo = 16384;

    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char out[CHUNKdo];

    vector<unsigned char> result;

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        throw exc :: error_from_deflateInit_t{};

    /* input data specified in one big chunk */
    strm.avail_in = src.size();
    flush = Z_FINISH; // instead of Z_NO_FLUSH, as this is the entire data
    strm.next_in = src.data();

    /* run deflate() on input until output buffer not full, finish
       compression if all of source has been read in */
    do {
        strm.avail_out = CHUNKdo;
        strm.next_out = out;
        ret = deflate(&strm, flush);    /* no bad return value */
        assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
        have = CHUNKdo - strm.avail_out;
        result.insert(result.end(), out, out+have);
    } while (strm.avail_out == 0);
    assert(strm.avail_in == 0);     /* all input will be used */

    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    cerr << "result.size() = " << result.size() << '\n';

    return result;
}

/* Decompress from file source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files. */
void inf(FILE *source, FILE *dest)
{
    constexpr int CHUNKii = 16384;
    constexpr int CHUNKio = 163;

    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNKii];
    unsigned char out[CHUNKio];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        throw exc:: error_from_inflateInit_t{ret};

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = fread(in, 1, CHUNKii, source);
        if (ferror(source)) {
            (void)inflateEnd(&strm);
            throw exc:: error_from_inflateInit_t{Z_ERRNO};
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNKio;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                throw exc:: error_from_inflateInit_t{ret};
            }
            have = CHUNKio - strm.avail_out;
            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
                (void)inflateEnd(&strm);
                throw exc:: error_from_inflateInit_t{Z_ERRNO};
            }
        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    if(ret != Z_STREAM_END) {
        throw exc:: error_from_inflateInit_t{Z_DATA_ERROR};
    }
}

/* report a zlib or i/o error */
void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
    case Z_ERRNO:
        if (ferror(stdin))
            fputs("error reading stdin\n", stderr);
        if (ferror(stdout))
            fputs("error writing stdout\n", stderr);
        break;
    case Z_STREAM_ERROR:
        fputs("invalid compression level\n", stderr);
        break;
    case Z_DATA_ERROR:
        fputs("invalid or incomplete deflate data\n", stderr);
        break;
    case Z_MEM_ERROR:
        fputs("out of memory\n", stderr);
        break;
    case Z_VERSION_ERROR:
        fputs("zlib version mismatch!\n", stderr);
    }
}

/* compress or decompress from stdin to stdout */
int main(int argc, char **argv)
{
    /* avoid end-of-line conversions */
    SET_BINARY_MODE(stdin);
    SET_BINARY_MODE(stdout);

    /* do compression if no arguments */
    if (argc == 1) {
        auto input_data = read_FILE_into_vector(stdin);
        vector<unsigned char> result = def(std::move(input_data), Z_DEFAULT_COMPRESSION);
        std:: cout.write    (   reinterpret_cast<char*>(result.data())
                            ,   result.size()
                            );
        return 0;
    }

    /* do decompression if -d specified */
    else if (argc == 2 && strcmp(argv[1], "-d") == 0) {
        inf(stdin, stdout);
        return 0;
    }

    /* otherwise, report usage */
    else {
        fputs("zpipe usage: zpipe [-d] < source > dest\n", stderr);
        return 1;
    }
}
