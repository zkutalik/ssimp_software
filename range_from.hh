#ifndef AMD_RANGE_FROM_HH
#define AMD_RANGE_FROM_HH

#include"range.hh"
namespace range {
namespace from {

template<typename T_vl>
struct from_ifstream_t {
    T_vl    m_f;
    std:: string  m_current;

    static_assert(!std:: is_rvalue_reference<T_vl>{} ,"");
    using value_type = std:: string;

    template<typename T>
    from_ifstream_t(T&& f)      :   m_f(AMD_FORWARD(f))
                                ,   m_current{}
    {
        advance();
    }

    std:: string            front_val()     const   { return m_current; }
    void                    advance()               {
        std:: getline(m_f, m_current);
    }
    bool                    empty()         const   { return false; }
};

template<typename T_vl>
auto    ifstream(T_vl && f) {
    return from_ifstream_t<T_vl>( AMD_FORWARD(f) );
}

} //namespace from
} //namespace range

#endif
