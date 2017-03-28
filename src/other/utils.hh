#include <cstdint>
#include <ostream>
#include <vector>

namespace utils {

template <typename T>
std:: intmax_t ssize(T const & x) {
    return x.size();
}

template<typename T>
std::ostream & operator<< (std:: ostream &o, const std:: vector<T> &v) {
    o << '[';
    for(auto & e : v) {
        if(&e != &v.front())
            o << ',';
        o << e;
    }
    o << ']';
    return o;
}

} // namespace utils
