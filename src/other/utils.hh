#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <utility> // for tuple::get

namespace utils {
namespace impl {
    template<typename ...T>
    struct voider { using type = void; };
}

template<typename ...T>
using void_t = typename impl:: voider<T...>:: type;

namespace impl {
    template<typename, template<class...> class Template, typename ...Args>
    struct can_apply_impl
        : public std:: false_type {};

    template<template<class...> class Template, typename ...Args>
    struct can_apply_impl<void_t<Template<Args...>>, Template, Args...>
        : public std:: true_type {
            using res = Template<Args...>;
    };
}

template<template<class...> class Template, typename ...Args>
using can_apply = impl:: can_apply_impl<void, Template, Args...>;

template<typename F, typename G, typename H>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F, G, H> &pr) {
    o << '('
        << std::get<0>(pr)
        << ','
        << std::get<1>(pr)
        << ','
        << std::get<2>(pr)
        << ')';
    return o;
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
template<typename T, typename ...Ts>
auto mk_vector(T t, Ts ...ts) -> std:: vector<T> {
    return std:: vector<T>{t, ts...};
}
template<typename T>
void print_type(T&&) {
    std:: cout << __PRETTY_FUNCTION__ << std:: endl;
}

template<typename T>
void pop_back_expected(std:: vector<T> &v, T const & expected) {
    assert(!v.empty());
    assert(expected == v.back());
    v.pop_back();
}
inline
void pop_back_expected_or_nan(std:: vector<long double> &v, long double expected) {
    assert(!v.empty());
    if(std::isnan(v.back())) {
        v.pop_back();
        return;
    }
    else
        pop_back_expected(v, expected);
}
template<typename T>
T    pop_back_and_return(std:: vector<T> &v) {
    assert(!v.empty());
    T was_at_the_back = std::move(v.back());
    v.pop_back();
    return was_at_the_back;
}
template<typename T>
std::ptrdiff_t ssize(const T& t) {
    return t.size();
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
std:: vector<std:: string>   tokenize       ( std:: string      const & line
                                            , char                delimiter
        );

} // namespace utils
