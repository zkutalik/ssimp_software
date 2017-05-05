/* This file is based on 'zpipe.c' from Mark Adler.
 * Here is a little on that file:
 *       zpipe.c: example of proper use of zlib's inflate() and deflate()
 *       Not copyrighted -- provided to the public domain
 *       Version 1.4  11 December 2005  Mark Adler
 * */

#include"zlib-vector.hh"
#include"../bits.and.pieces/ASSERT.hh"

using zlib_vector:: vec_t;

#include<zlib.h>
#include<cassert>

namespace exc {
    // Too lazy to make proper std::exception objects to throw. So
    // I'll just throw these instead:
    struct error_from_deflateInit_t {};
    struct error_from_inflateInit_t {
        decltype(Z_DATA_ERROR) m_error_code;
    };
}

vec_t zlib_vector:: deflate(vec_t src)
{
    int level = 9; // for best compression
    vec_t result;

    constexpr int CHUNK = 16384;

    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char out[CHUNK];

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
        strm.avail_out = CHUNK;
        strm.next_out = out;
        ret = deflate(&strm, flush);    /* no bad return value */
        assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
        have = CHUNK - strm.avail_out;
        result.insert(result.end(), out, out+have);
    } while (strm.avail_out == 0);
    assert(strm.avail_in == 0);     /* all input will be used */

    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);

    return result;
}
vec_t zlib_vector:: inflate(vec_t src)
{
    vec_t result;
    constexpr int CHUNK = 16384;

    int ret;
    unsigned have;
    z_stream strm;
    unsigned char out[CHUNK];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        throw exc:: error_from_inflateInit_t{ret};

    strm.avail_in = src.size();
    strm.next_in = src.data();

    /* run inflate() on input until output buffer not full */
    do {
        strm.avail_out = CHUNK;
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
        have = CHUNK - strm.avail_out;
        result.insert(result.end(), out, out+have);
    } while (strm.avail_out == 0);

    // We deflated one big chunk of input, so ret
    // will be Z_STREAM_END now
    assert(ret == Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return result;
}
