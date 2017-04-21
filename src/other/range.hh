#pragma once
#include<type_traits>
#include<utility>
#include<stdexcept>
#include<numeric> // for std::accumulate

#include<vector>
#include<tuple>
#include<algorithm>

#include<memory>
#include<iostream> // just so I can implement operator <<

#include"utils.hh"

#define AMD_RANGE_DECLTYPE_AND_RETURN( expr )    decltype( expr ) { return expr ; }
#define AMD_FORWARD(                   expr )    std::forward<decltype(expr)>(expr)

using utils:: can_apply;
using utils:: void_t;

namespace range {
    struct range_tag {}; // *all* ranges will inherit this, if nothing else
    template<typename R>
    using is_of_range_tag = std:: enable_if_t<std::is_base_of<range_tag, std::remove_reference_t<R>>{}>;
    // A range is either non-owning:
    //   - initialized with lvalue container (or an lvalue range)
    //   - the range itself is copyable
    // or, is owning:
    //   - move-constructs it's argument in (to a unique_ptr perhaps)
    //   - the range cannot be (implicitly) copied

    //   the value_type may not be a reference

    /* Possible methods:
     *
     *   bool empty()      - no more input can be read, or output written
     *   T    pull()       - like from an input stream (front_val followed by advance) {might throw pull_from_empty_range_error}
     *   T    front_val()  - repeated reads read from the same location
     *   T&   front_ref()  - repeated reads read from the same location
     *
     *   void advance()    - skip the current value
     */

    // First, I'll give an example of a simple range-like object. What could be simpler than a pair of iterators?
    template<typename b_t, typename e_t>
    struct range_from_begin_end_t : public range_tag {
        b_t m_b;
        e_t m_e;

        using value_type = std::decay_t<decltype(*m_b)>;

        range_from_begin_end_t(b_t b, e_t e) : m_b(move(b)), m_e(move(e)) {}

        auto begin() const { return m_b; }
        auto end  () const { return m_e; }
        bool empty() const { return m_b == m_e; }
        void advance()     { ++m_b; }
        auto current_it() const { return m_b; }
        // should consider a more flexible front_ref that tries to
        // return the least cv-qualified version that it can
        value_type      const &  front_ref() const   { return *m_b; }
    };
    template<typename b_t, typename e_t>
    range_from_begin_end_t<b_t,e_t> range_from_begin_end(b_t b, e_t e) {
        return {move(b), move(e)};
    }
    template<typename C>
    auto range_from_begin_end(C &c) -> AMD_RANGE_DECLTYPE_AND_RETURN( range_from_begin_end(c.begin(), c.end()) )

    template<typename I>
    struct range_ints_t : public range_tag {
        struct iter_is_own_value {
            I m_i;

            bool    operator!=  (iter_is_own_value const & other) const {
                return m_i != other.m_i;
            }
            void    operator++  ()          { ++m_i; }
            I     operator*   () const    { return m_i; }
        };
        I m_b;
        I m_e;

        using value_type = I;

        range_ints_t(I e)        : m_b(0), m_e(e) {}
        range_ints_t(I b, I e) : m_b(b), m_e(e) {}

        iter_is_own_value begin() const { return {m_b}; }
        iter_is_own_value end  () const { return {m_e}; }
        bool        empty()     const   { return  m_b == m_e; }
        void        advance()           {       ++m_b; }
        I           front_val  () const { return  m_b; }
    };
    inline
    range_ints_t<int> ints(int e) {
        return {e};
    }

    struct pull_from_empty_range_error : public std:: runtime_error {
        pull_from_empty_range_error() : std:: runtime_error("attempted pull() from an empty range") {}
    };

    template<typename V>
    struct from_vector_t
        // If V is *not* a reference, we want this to be non-copyable
        : public std:: conditional_t<  std:: is_reference<V>{}, utils:: empty, utils:: non_copyable_empty >
    { // V may, not may not, be a reference
        V m_v;
        size_t m_i;

        static_assert( !std:: is_rvalue_reference<V>{} ,"");

        template<typename V2>
        from_vector_t(V2 &&v, size_t i) : m_v(AMD_FORWARD(v)), m_i(i) {}

        bool    empty()         const   { return m_i >= m_v.size(); }
        void    advance  ()             { ++m_i; }
        decltype(auto)    front_ref()     const   { return get_fwd().at(m_i); }
        decltype(auto)    front_ref()             { return get_fwd().at(m_i); }


        decltype(auto) get_fwd() {
            return std::forward<V>(m_v);
        }
        decltype(auto) get_fwd() const {
            // We need to add const to V before forwarding
            // I wish std::forward didn't make things which
            // were 'less const'.

            // Anyway, V is either non-ref or l-ref, but not r-ref
            static_assert(!std:: is_rvalue_reference<V>{} ,"");

            using V_const = utils:: add_const_even_if_ref_t<V>;
            using fwd_t_const = decltype(std::forward< V_const >(m_v));
            static_assert( std:: is_same<V_const &&, fwd_t_const> {} ,"");

            // Check V and V_const have the same 'referencedness'
            static_assert( std:: is_lvalue_reference<V_const>{} == std:: is_lvalue_reference<V> {} ,"");
            static_assert( std:: is_rvalue_reference<V_const>{} == std:: is_rvalue_reference<V> {} ,"");

            // ... but different 'constness'
            static_assert(!std:: is_const<std::remove_reference_t<V>> {} ,"");
            static_assert( std:: is_const<std::remove_reference_t<V_const>> {} ,"");
            static_assert(!std:: is_same<V_const, V> {} ,"");
            static_assert( std:: is_same<   std::decay_t<V_const>
                                        ,   std::decay_t<V>
                                        > {} ,"");
            return std::forward<V_const>(m_v);
        }

        template<typename T>
        auto push_back(T &&t) -> AMD_RANGE_DECLTYPE_AND_RETURN( get_fwd().push_back( AMD_FORWARD(t) ) );

        template<typename T>
        auto push_back(T &&t) const -> AMD_RANGE_DECLTYPE_AND_RETURN( get_fwd().push_back( AMD_FORWARD(t) ) );
    };

    template<typename V>
    auto from_vector(V &&v) {
        return from_vector_t<V>
                ( std::forward<V>(v), 0 );
    }

    template<typename R>
    auto front_ref(R&& r)
    -> decltype( std::forward<R>(r).front_ref() )
    {
        return   std::forward<R>(r).front_ref();
    }

    template<typename R, typename T>
    auto push_back(R&& r, T &&t)
    -> decltype( std::forward<R>(r).push_back(std::forward<T>(t)) )
    {
        return   std::forward<R>(r).push_back(std::forward<T>(t));
    }

    template<typename R >
    auto front_val_impl(R&& r, utils:: priority_tag<5>) -> AMD_RANGE_DECLTYPE_AND_RETURN(
            std::forward<R>(r).front_val() )
    template<typename R >
    auto front_val_impl(R&& r, utils:: priority_tag<3>) -> AMD_RANGE_DECLTYPE_AND_RETURN(
            utils:: un_lref( std::forward<R>(r).front_ref() ) )
    template<typename R >
    auto front_val_impl(R&& r, utils:: priority_tag<1>)
    -> std:: decay_t<decltype( *std::forward<R>(r).begin() )>
    {
        return *std::forward<R>(r).begin();
    }

    // Next, synthesize range::empty, for those with an empty() method
    template<typename R >
    auto empty(R&& r) -> AMD_RANGE_DECLTYPE_AND_RETURN( std::forward<R>(r).empty() )

    // Synthesize range:: advance
    template<typename R >
    auto advance(R&& r) -> AMD_RANGE_DECLTYPE_AND_RETURN( std::forward<R>(r).advance() )

    template<typename idx, typename ... range_types>
    struct zip_val_t;
    template<size_t ...Is, typename ...Rs>
    struct zip_val_t<std::index_sequence<Is...>, Rs...> : public range_tag {
        std:: tuple<Rs...> m_ranges;

        using value_type = std:: tuple< decltype( front_val(std::get<Is>(m_ranges)) ) ... >;

        static_assert( utils:: and_all( std:: is_same<Rs, std::decay_t<Rs> >{}... ), "");

        template<typename ...Ts>
        zip_val_t(Ts&& ...ranges)
            : m_ranges( std::forward<decltype(ranges)>(ranges)... )
            {}

        value_type          front_val()     const;

        bool                empty()         const {
            bool e0 = range:: empty(std::get<0>(m_ranges));
            bool all_the_same = utils:: and_all( e0 == range:: empty(std::get<Is>(m_ranges)) ...);
            assert(all_the_same);
            return e0;
        }
        void                advance() {
            auto ignore_me = {((void)std::get<Is>(m_ranges).advance(),false)...};
            (void)ignore_me;
        }

        template<typename id = utils:: id>
        auto push_back(value_type p)
        -> AMD_RANGE_DECLTYPE_AND_RETURN
        (
                   utils:: ignore( ((void)range:: push_back(std::get<Is>( id{}(m_ranges)), std::get<Is>(p)),0) ... )
        )
    };
    template<typename R >
    auto front_val(R&& r) -> AMD_RANGE_DECLTYPE_AND_RETURN( front_val_impl(std::forward<R>(r), utils:: priority_tag<9>()) );
    template<size_t ...Is, typename ...Rs>
    typename zip_val_t<std::index_sequence<Is...>, Rs...> ::
        value_type
    zip_val_t<std::index_sequence<Is...>, Rs...> ::
        front_val()     const {
            return  value_type  { range:: front_val(std::get<Is>(m_ranges)) ...  };
        }

    template<typename ...Rs>
    zip_val_t< std:: index_sequence_for<Rs...> ,std::decay_t<Rs>... >
    zip_val(Rs && ...ranges) {
        return { std::forward<decltype(ranges)>(ranges)... };
    }

    namespace impl {
    template<typename R>
    auto pull_impl(R&& r, utils:: priority_tag<3>) -> AMD_RANGE_DECLTYPE_AND_RETURN (
        std:: forward<R>(r).pull()
    )
    template<typename R>
    auto pull_impl(R&& r, utils:: priority_tag<2>) -> decltype( range:: front_val( std:: forward<R>(r) ) )
    {
        auto x = range:: front_val( std:: forward<R>(r) );
                 range:: advance  ( std:: forward<R>(r) );
        return x;
    }
    }

    template<typename R>
    auto pull(R&& r) {
        return impl:: pull_impl(std::forward<R>(r), utils:: priority_tag<9>{});
    }

    template<typename R>
    struct begin_end_for_range_for {

        R * m_range_pointer;

        bool operator != ( const begin_end_for_range_for & other ) const {
            if(other.m_range_pointer == nullptr && this->m_range_pointer == nullptr) {
                // both are end-iterators.
                return false;
            }
            if(other.m_range_pointer == nullptr && this->m_range_pointer != nullptr) {
                return !this->m_range_pointer->empty();
            }
            if(other.m_range_pointer != nullptr && this->m_range_pointer == nullptr) {
                return !other.m_range_pointer->empty();
            }
            // BUG INCOMPLETE
            throw std:: runtime_error("range:: range-based-for (incomplete checks here)");
        }
        void error_now_if_I_am_empty() const {
              this->m_range_pointer           || [](){throw std:: runtime_error("range:: range-based-for (++end)"    );return false;}();
            (!this->m_range_pointer->empty()) || [](){throw std:: runtime_error("range:: range-based-for (++empty())");return false;}();
        }
        begin_end_for_range_for& operator++() {
            error_now_if_I_am_empty();
            this->m_range_pointer->advance();
            return *this;
        }

        template<class...  , typename Myself =              const R& >
        auto operator_star_impl(utils:: priority_tag<2>)    const
        -> decltype( std::declval<Myself>().front_ref() )
        {
            error_now_if_I_am_empty();
            return   this->m_range_pointer->front_ref();
        }

        template<class...  , typename Myself =              const R& >
        auto operator_star_impl(utils:: priority_tag<1>)    const
        -> decltype( std::declval<Myself>().front_val() )
        {
            error_now_if_I_am_empty();
            return   this->m_range_pointer->front_val();
        }
        auto operator*() const
            -> AMD_RANGE_DECLTYPE_AND_RETURN(
                   operator_star_impl(utils:: priority_tag<9>{})
            )
    };

    template<typename R>
    auto end  (R &)
    -> AMD_RANGE_DECLTYPE_AND_RETURN(begin_end_for_range_for<R> { nullptr })
    template<typename R >
    auto begin(R &r)
    -> AMD_RANGE_DECLTYPE_AND_RETURN(begin_end_for_range_for<R> { &r })

    // These next two are deliberately to cause linker failures.
    // Instead I should do the priority_tag thing to directly use
    // .begin and .end methods when present
    template<typename R>
    auto end  (R && r) -> decltype(r.end()  );
    template<typename R>
    auto begin(R && r) -> decltype(r.begin());

    template<typename R, class..., typename Rnonref = std::remove_reference_t<R>, typename = std:: enable_if_t<std::is_base_of<range_tag, Rnonref >{}> >
    decltype(auto) operator<< (std:: ostream &o, R r) {
        if(r.empty()) {
            o << "[]";
            return o;
        }
        o << '[';
        for(;;) {
            o << front_val(r);
            r.advance();
            if(r.empty())
                break;
            o << ',';
        }
        o << ']';
        return o;
    }

    template        <typename Rb, typename Re, typename G >
    auto sort(range_from_begin_end_t<Rb,Re> r, G &&init) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::sort( r.m_b , r.m_e , std::forward<G>(init)))
    template        <typename Rb, typename Re>
    auto max_element(range_from_begin_end_t<Rb,Re> r) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::max_element( r.m_b , r.m_e))
    template        <typename Rb, typename Re, typename G >
    auto shuffle(range_from_begin_end_t<Rb,Re> r, G &&g) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::shuffle( r.m_b , r.m_e , std::forward<G>(g)))
    template        <typename Rb, typename Re, typename G >
    auto accumulate(range_from_begin_end_t<Rb,Re> r, G &&init) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::accumulate( r.m_b , r.m_e , std::forward<G>(init)))
    template        <typename Rb, typename Re>
    auto next_permutation(range_from_begin_end_t<Rb,Re> r) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::next_permutation( r.m_b , r.m_e))

} // namespace range
