#ifndef AMD_UTILS_HH__
#define AMD_UTILS_HH__

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <utility> // for tuple::get
#include <vector>
#include <string>

#include "ASSERT.hh"

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
template<typename F>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F> &pr);
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
inline
auto nice_operator_shift_left(std:: string t)
-> std:: string
{
    std:: vector<char> q;
    q.reserve(2+2*t.size());
    q.push_back('"');
    for(char c : t) {
        switch(c) {
            break; case '"': q.push_back('\\'); q.push_back('"');
            break; case '\n': q.push_back('\\'); q.push_back('n');
            break; case '\t': q.push_back('\\'); q.push_back('t');
            break; case '\\': q.push_back('\\'); q.push_back('\\');
            break; default:
                q.push_back(c);
        }
    }
    q.push_back('"');
    return std::string{ q.begin(), q.end() };
}

template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: pair<F, G> &pr) {
    o << '('
        << nice_operator_shift_left( std::get<0>(pr) )
        << ','
        << nice_operator_shift_left( std::get<1>(pr) )
        << ')';
    return o;
}
template<typename F>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F> &pr) {
    o << '('
        << nice_operator_shift_left(std::get<0>(pr))
        << ')';
    return o;
}
template<typename F, typename G>
std::ostream & operator<< (std:: ostream &o, const std:: tuple<F, G> &pr) {
    o << '('
        << nice_operator_shift_left(std::get<0>(pr))
        << ','
        << nice_operator_shift_left(std::get<1>(pr))
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
template<typename ...Ts>
void print_type(Ts && ...) {
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
constexpr
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

// and_all
//auto and_all(void) = delete; // This seems to break gcc. Should consider a bug report
template<typename T>
constexpr
auto and_all(T t) { return t; }
template<typename T0, typename T1, typename ...Ts>
constexpr
auto and_all(T0 t0, T1 t1, Ts ...ts) { return and_all( t0&&t1, ts... ); }

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

struct save_ostream_briefly {
    std::ostream               & m_o;
    std::ios::fmtflags      m_f;

    save_ostream_briefly(std::ostream & o) : m_o(o) {
        m_f = m_o.flags();
    }
    ~save_ostream_briefly() {
        m_o.flags(m_f);
    }
};


//   make_a_pack_and_apply_it(F&&f);
namespace detail {
    template<typename T, typename F, T ...Idxs>
    decltype(auto) make_a_pack_and_apply_it(std:: integer_sequence<T, Idxs...>, F&&f) {
        return std::forward<F>(f) ( std:: integral_constant<T, Idxs>{} ... );
    }
}
template<size_t N, typename T, typename F>
decltype(auto) make_a_pack_and_apply_it(F&&f) {
    return  detail:: make_a_pack_and_apply_it<T>
            (   std:: make_index_sequence<N>{}
            ,   std::forward<F>(f)
            );
}

template<char ... chars>
struct char_pack {
    constexpr static char   c_str0_[] =  {chars..., '\0'};

    constexpr static size_t size()      { return sizeof...(chars); }
    constexpr static char   const   (&c_str0(void)) [size()+1]     { return c_str0_; }
    constexpr static char   at(size_t i)        { return c_str0_[i]; }
    template<typename ...C>
    constexpr static size_t find_first_of(C ... targets) {
        for(size_t i=0; i<   1+ size(); ++i) {
            if( !utils:: and_all( at(i) != targets ... ))
                return i;
        }
        return -1;
    }
    template<size_t l>
    constexpr static auto   substr(void) {
        return
        make_a_pack_and_apply_it<l, size_t>([](auto ... idxs) {
            return utils:: char_pack< char_pack:: at(idxs) ... >{};
        });
    }
    template<size_t b, size_t e>
    constexpr static auto   substr(void) {
        static_assert( e>=b ,"");
        return
        make_a_pack_and_apply_it<e-b, size_t>([](auto ... idxs) {
            return utils:: char_pack< char_pack:: at(b+idxs) ... >{};
        });
    }
};
template<char ... chars>
constexpr char   char_pack<chars...>:: c_str0_[];

template<typename T, T c>
struct compile_time_constant_as_a_type {
    constexpr
    bool    operator==  (T other)   const   {   return c == other; }

    constexpr
    operator T  ()                  const   {   return c; }
};

template<typename T, T c>
constexpr compile_time_constant_as_a_type<T, c>     cx_val  = {}; // Consider using my own type here instead of integral_constant?

template<typename T, T c1, T c2>
constexpr
auto    operator+   (   compile_time_constant_as_a_type<T,c1>
                    ,   compile_time_constant_as_a_type<T,c2>   )
{
    return cx_val<T, c1+c2>;
}
template<typename T, T c1, T c2>
constexpr
auto    operator==  (   compile_time_constant_as_a_type<T,c1>
                    ,   compile_time_constant_as_a_type<T,c2>   )
{
    return cx_val<bool, c1==c2>;
}

template<typename T, typename U>
bool is_in(T && t, std:: initializer_list<U> l) {
    for(auto && x : l) {
        if(t==x)
            return true;
    }
    return false;
}

inline
bool startsWith(std:: string const & s, std:: string const & prefix) {
    return s.substr(0, prefix.length()) == prefix;
}

} // namespace utils

#endif
