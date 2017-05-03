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

/* Next few lines are needed to help ADL - they forward *as-is* */
template<typename T> auto  begin        (T&&t)  -> AMD_RANGE_DECLTYPE_AND_RETURN( range:: begin     (AMD_FORWARD(t)) )
template<typename T> auto  end          (T&&t)  -> AMD_RANGE_DECLTYPE_AND_RETURN( range:: end       (AMD_FORWARD(t)) )
template<typename T> auto  front_val    (T&&t)  -> AMD_RANGE_DECLTYPE_AND_RETURN( range:: front_val (AMD_FORWARD(t)) )

} //namespace from
} //namespace range

#endif
