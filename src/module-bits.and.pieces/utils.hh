#ifndef AMD_UTILS_HH__
#define AMD_UTILS_HH__

#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <utility> // for tuple::get
#include <vector>
#include <string>

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

template<typename T>
decltype(auto) consider_quoting(T &&t) {
    return std::forward<T>(t);
}
//template<>
inline
std:: string consider_quoting(std::string const &t) {
    std:: ostringstream oss;
    oss << '"';
    for(char c : t) {
        switch(c) {
            break; case '\t': oss << "\\t";
            break; case '"' : oss << "\\\"";
            break; case '\\' : oss << "\\\\";
            break; default  : oss << c;
        }
    }
    oss << '"';
    return oss.str();
}

template<typename T>
std::ostream & operator<< (std:: ostream &o, const std:: vector<T> &v);
template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: pair<F, G> &pr);
template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F, G> &pr);
template<typename F, typename G, typename H>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F, G, H> &pr);
template<typename T>
auto operator<< (std:: ostream &o, T &t)
-> decltype(  o << t.to_string() );

template<typename T>
auto nice_operator_shift_left(T &&t)
-> decltype( std:: forward<T>(t)){ return std:: forward<T>(t); }
inline
auto nice_operator_shift_left(uint8_t t)
-> int
{ return t; }
inline
auto nice_operator_shift_left(bool t)
-> char
{ return t?'T':'F'; }

template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: pair<F, G> &pr) {
    o << '('
        << nice_operator_shift_left( std::get<0>(pr) )
        << ','
        << nice_operator_shift_left( std::get<1>(pr) )
        << ')';
    return o;
}
template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F, G> &pr) {
    o << '('
        << std::get<0>(pr)
        << ','
        << std::get<1>(pr)
        << ')';
    return o;
}
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
template<>
inline
std::ostream & operator<< <bool> (std:: ostream &o, const std:: vector<bool> &v) {
    o << '[';
    for(size_t i = 0; i<v.size(); ++i) {
        if(i!=0)
            o << ',';
        o << nice_operator_shift_left(v.at(i));
    }
    o << ']';


    return o;
}
template<typename T>
auto operator<< (std:: ostream &o, T &t)
-> decltype(  o << t.to_string() )
{
    return o << t.to_string();
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
void print_type() {
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

inline
std:: vector<std:: string>   tokenize       ( std:: string      const & line
                                            , char                      delimiter
        ) {
    int pos = 0;
    std:: vector<std:: string> fields;
    while(1) {
        auto next_delim = line.find(delimiter, pos);
        fields.push_back( line.substr(pos, next_delim-pos) );
        if(next_delim== std:: string:: npos)
            break;
        pos=next_delim+1;
    }
    return fields;
}

template<int i>
struct priority_tag;
template<int i>
struct priority_tag : public priority_tag<i-1> {};
template<>
struct priority_tag<0> {};

template<typename ...>
struct and_all_impl;

template<typename ...Ts>
constexpr
auto and_all(Ts ...ts) { return and_all_impl<Ts...>::impl(ts...); }

template<typename T0>
struct and_all_impl<T0> {
    static constexpr
    T0  impl(T0 t0) { return t0; }
};
template<typename T0, typename T1, typename ...Ts>
struct and_all_impl<T0, T1, Ts...> {
    static constexpr
    auto  impl(T0 t0, T1 t1, Ts ...ts) { return and_all(t0 && t1, ts...); }
};

// This is my own 'apply', as it's officially only in c++17, not c++14
namespace detail {
template <class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F &&f, Tuple &&t, std::index_sequence<I...>)
{
    return std::forward<F>(f)(std::get<I>(std::forward<Tuple>(t))...);
}
}  // namespace detail

template <class F, class Tuple>
constexpr decltype(auto) apply(F &&f, Tuple &&t)
{
    return detail::apply_impl(
        std::forward<F>(f), std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{}>{});
}

struct non_copyable_empty {
    non_copyable_empty             ()                           = default;
    non_copyable_empty             (non_copyable_empty const &) = delete;
    non_copyable_empty             (non_copyable_empty      &&) = default;
    non_copyable_empty & operator= (non_copyable_empty const &) = delete;
    non_copyable_empty & operator= (non_copyable_empty      &&) = default;
};
struct empty {
};

template<typename ...Ts>
void ignore(Ts && ...) {}

struct id {
    template<typename T>
    decltype(auto) operator() (T &&t) { return std::forward<T>(t); }
};

template<typename T>
struct add_const_under_ref {
    static_assert( !std:: is_reference<T>{} ,"");
    // don't define 'type' here, as this refuses
    // to add const to a non-ref type
};
template<typename T>
struct add_const_under_ref<T&> {
    static_assert( !std:: is_reference<T>{} ,"");
    using type = std:: add_const_t<T> &;
};
template<typename T>
struct add_const_under_ref<T&&> {
    static_assert( !std:: is_reference<T>{} ,"");
    using type = std:: add_const_t<T> &&;
};
template<typename T>
using add_const_under_ref_t = typename add_const_under_ref<T>::type;

template<typename T>
struct add_const_even_if_ref {
    static_assert( !std:: is_reference<T>{} ,"");
    using type = std:: add_const_t<T>;
};
template<typename T>
struct add_const_even_if_ref<T&> {
    static_assert( !std:: is_reference<T>{} ,"");
    using type = std:: add_const_t<T> &;
};
template<typename T>
struct add_const_even_if_ref<T&&> {
    static_assert( !std:: is_reference<T>{} ,"");
    using type = std:: add_const_t<T> &&;
};
template<typename T>
using add_const_even_if_ref_t = typename add_const_even_if_ref<T>::type;

template<typename T>
T   un_lref(T &t) { return t; }
template<typename T>
T   un_lref(T &&) = delete;

inline double ELAPSED(void) {
       return double(clock()) / CLOCKS_PER_SEC;
}

struct {} stdget0;
struct {} stdget1;
template<typename T>
decltype(auto) operator| (T&& t, decltype(stdget0)) {
    return std:: get<0>( std::forward<T>(t) );
}
template<typename T>
decltype(auto) operator| (T&& t, decltype(stdget1)) {
    return std:: get<1>( std::forward<T>(t) );
}


} // namespace utils

#endif
