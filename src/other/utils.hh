#include <cstdint>
#include <ostream>
#include <iostream>
#include <stdexcept>
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
template<>
inline
std::ostream & operator<< <uint8_t> (std:: ostream &o, const std:: vector<uint8_t> &v) {
    o << '[';
    for(auto & e : v) {
        if(&e != &v.front())
            o << ',';
        o << (int)e;
    }
    o << ']';
    return o;
}

template<typename T>
T lexical_cast(std:: string const & s);

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
template<>
inline
double lexical_cast<double>(std:: string const & s) {
    double val;
    int n;
    int ret = sscanf(s.c_str(), "%lg %n", &val, &n);
    if(ret == 1 && n == ssize(s)) {
        return val;
    }
    else
        throw std:: invalid_argument{std::string("Can't parse this double [") + s + ']'};
}

template<typename T>
void print_type(T&&) {
    std:: cout << __PRETTY_FUNCTION__ << std:: endl;
}

} // namespace utils
