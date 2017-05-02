#include<vector>
namespace zlib_vector {
    using vec_t = std:: vector<unsigned char>;

    vec_t deflate(vec_t src);
    vec_t inflate(vec_t src);

}

