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

template<typename T>
int lexical_cast(std:: string const & s);

template<>
inline
int lexical_cast<int>(std:: string const & s) {
    int d;
    int n;
    int ret = sscanf(s.c_str(), "%d %n", &d, &n);
    if(ret == 1 && n == ssize(s)) {
        return d;
    }
    else
        throw std:: invalid_argument{"can't parse this int"};
}

} // namespace utils
