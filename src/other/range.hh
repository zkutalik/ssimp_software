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
//#include"little.hh"
#include "void_t_and_related.hh"

#define AMD_RANGE_DECLTYPE_AND_RETURN( expr )    decltype( expr ) { return expr ; }

using void_t_and_related:: can_apply;
using void_t_and_related:: void_t;

namespace range {
    struct range_tag {}; // *all* ranges will inherit this, if nothing else
    template<typename R>
    using is_of_range_tag = std:: enable_if_t<std::is_base_of<range_tag, std::remove_reference_t<R>>{}>;
    // A range is either non-owning:
    //   - initialized with lvalue container (or an lvalue range)
    //   - the range itself is non-copyable
    // or, is owning:
    //   - move-constructs it's argument in
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
    struct pull_from_empty_range_error : public std:: runtime_error {
        pull_from_empty_range_error() : std:: runtime_error("attempted pull() from an empty range") {}
    };

    // test whether the object has particular methods
    namespace {
        template<typename R> using tester_lval_has_method_front_ref = decltype( std:: declval<R&>().front_ref() );
        template<typename R> using tester_lval_has_method_front_val = decltype( std:: declval<R&>().front_val() );
        template<typename R> using tester_rval_has_method_pull      = decltype( std:: declval<R>().pull() );
        template<typename R> using tester_has_front_ref        = decltype( front_ref(std:: declval<R&>()) );
        template<typename R> using tester_has_front_val        = decltype( front_val(std:: declval<R&>()) );
        template<typename R> using tester_has_method_is_infinite  = decltype( std:: declval<R&>().is_infinite() );
    }
    // in the next few tests, I should think more about lval/rval of the range itself
    template<typename R> using has_method_front_ref = can_apply<tester_lval_has_method_front_ref,R>;
    template<typename R> using has_method_front_val = can_apply<tester_lval_has_method_front_val,R>;
    template<typename R> using has_method_pull      = can_apply<tester_rval_has_method_pull,R>;
    template<typename R> using has_front_ref        = can_apply<tester_has_front_ref,R>;
    template<typename R> using has_front_val        = can_apply<tester_has_front_val,R>;
    template<typename R> using has_method_is_infinite = can_apply<tester_has_method_is_infinite,R>;

    /* Now a few free functions, for use with zip at first */
    template<typename R
            , class ...
            , typename = std::enable_if_t< has_method_front_ref<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto front_ref(R&& r)
        -> decltype( std::forward<R>(r).front_ref() )
    {
        return std::forward<R>(r).front_ref();
    }
    template<typename R
            , class ...
            , typename = std::enable_if_t< has_method_front_val<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto front_val(R&& r)
        -> decltype( std::forward<R>(r).front_val() )
    {
        return std::forward<R>(r).front_val();
    }

    template<typename R
            , class ...
            , typename = std::enable_if_t< has_front_ref<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto front_ref_then_val(R&& r)
    -> decltype( front_ref(std::forward<R>(r)) )
    {
        return front_ref(std::forward<R>(r));
    }
    template<typename R
            , class ...
            , typename = std::enable_if_t< !has_front_ref<R>{} && has_front_val<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto front_ref_then_val(R&& r)
    -> decltype( front_val(std::forward<R>(r)) )
    {
        return front_val(std::forward<R>(r));
    }

    template<typename R
            , class ...
            , typename = std::enable_if_t< has_method_pull<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto pull(R&& r)
        -> decltype( std::forward<R>(r).pull() )
    {
              return std::forward<R>(r).pull();
    }
    template<typename R
            , class ...
            , typename = std::enable_if_t< !has_method_pull<R>{} >
            , typename = is_of_range_tag<R>
            >
    auto pull(R&& r)
        -> decltype( std::forward<R>(r).front_val() )
    {
              auto val = std::forward<R>(r).front_val();
              std::forward<R>(r).advance();
              return val;
    }

    template<typename R
            , class ...
            , typename = is_of_range_tag<R>
            , typename = std::enable_if_t< has_method_is_infinite<R>{} >
            , typename = void // just to disambig with the next one
            >
    constexpr bool is_infinite(R&& r) {
        return r.is_infinite();
    }
    template<typename R
            , class ...
            , typename = is_of_range_tag<R>
            , typename = std::enable_if_t<!has_method_is_infinite<R>{} >
            >
    constexpr bool is_infinite(R&&) {
        return false;
    }

    template<typename getter, typename underlying_type, typename underlying_refs_type, size_t ...is>
    struct zipped_range : public range_tag {
        private:
            underlying_type underlying;
        public:

            using value_type = std:: tuple<
                //typename std:: tuple_element<is, underlying_type>::type ::value_type & ...
                decltype(getter:: getter(std::get<is>(underlying))) ...
            >;

            zipped_range(underlying_refs_type &&underlying_) : underlying( std:: move(underlying_)) {
                // constructing a tuple of non-refs from a tuple of [lr]v-refs
            }

            bool empty() const {
                bool are_any_empty = std:: max( {std::get<is>(underlying).empty() ...} );
                if(! are_any_empty) {
                    // they're all non-empty, just return
                }
                else {
                    // at least one is empty. return true.
                    // But check first that they are *all* empty or infinite
                    bool all_are_empty_or_infinite = std:: min( {is_infinite(std::get<is>(underlying))||std::get<is>(underlying).empty() ... } );
                    all_are_empty_or_infinite || [](){throw std:: runtime_error("range:: imbalanced sources in zip range");return false;}();
                }
                return are_any_empty;
            }
            template<class ..., typename = void_t<decltype( getter:: getter( std::get<is>(underlying) ) ) ... >>
            value_type front_val() const {
                (!this->empty()) || [](){throw std:: runtime_error("range:: front_val() called on empty zip range");return false;}();
                return value_type{ getter:: getter(std::get<is>(underlying)) ... };
            }
            void advance() {
                int just_for_side_effects[] = { ((void)
                        std::get<is>(underlying).advance()
                        ,0)... };
                (void)just_for_side_effects;
            }

    };
    template<typename getter, size_t ...is, typename ...RangesRef>
    auto impl_zip(
            std:: integer_sequence<size_t, is...>
            ,std::tuple< RangesRef ... > &&underlying_refs
            ) {

        using underlying_refs_type = std:: remove_reference_t< decltype(underlying_refs) >;
        using underlying_type = std:: tuple<
            std:: remove_reference_t<
                typename std:: tuple_element<is, underlying_refs_type>::type
                > ... >;

        return zipped_range<getter, underlying_type, underlying_refs_type, is...>{ std:: move(underlying_refs) };
    }
    struct get_the_ref_not_val {
        template<typename R>
        static
        auto getter(R&& r)
            -> decltype( front_ref( std::forward<R>(r) ) )
        {
            return front_ref( std::forward<R>(r) );
        }
    };
    struct get_the_val_not_ref {
        template<typename R>
        static
        auto getter(R&& r)
            -> decltype( front_val( std::forward<R>(r) ) )
        {
            return front_val( std::forward<R>(r) );
        }
    };
    struct get_the_ref_OR_val {
        template<typename R>
        static
        auto getter(R&& r)
            -> decltype( front_ref_then_val( std::forward<R>(r) ) )
        {
            return front_ref_then_val( std::forward<R>(r) );
        }
    };
    template<typename ...Ranges
        , class..., typename = void_t< has_method_front_ref<Ranges> ...  >
        >
    auto zip_the_refs(Ranges&& ...ranges) {
        return impl_zip<get_the_ref_not_val>(
                std:: make_index_sequence< sizeof...(ranges) >{}
               ,std:: tuple< Ranges&&... >{ std::forward<Ranges>(ranges) ... }
               );
    }
    template<typename ...Ranges
        , class..., typename = void_t< has_method_front_val<Ranges> ...  >
        >
    auto zip_the_vals(Ranges&& ...ranges) {
        return impl_zip<get_the_val_not_ref>(
                std:: make_index_sequence< sizeof...(ranges) >{}
               ,std:: tuple< Ranges&&... >{ std::forward<Ranges>(ranges) ... }
               );
    }
    template<typename ...Ranges
        , class..., typename = void_t< has_method_front_val<Ranges> ...  >
        >
    auto zip(Ranges&& ...ranges)
    -> decltype(
               impl_zip<get_the_ref_OR_val>(
                std:: make_index_sequence< sizeof...(ranges) >{}
               ,std:: tuple< Ranges&&... >{ std::forward<Ranges>(ranges) ... }
               )
    )
    {
        return impl_zip<get_the_ref_OR_val>(
                std:: make_index_sequence< sizeof...(ranges) >{}
               ,std:: tuple< Ranges&&... >{ std::forward<Ranges>(ranges) ... }
               );
    }
    template<typename R>
    struct begin_end_for_range_for {
        std:: unique_ptr<R> m_range_pointer;
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
        begin_end_for_range_for& operator++() {
            this->m_range_pointer             || [](){throw std:: runtime_error("range:: range-based-for (++end)"    );return false;}();
            (!this->m_range_pointer->empty()) || [](){throw std:: runtime_error("range:: range-based-for (++empty())");return false;}();
            this->m_range_pointer->advance();
            return *this;
        }

        // Of the following two, only will is allowed to exist.
        // Same signature, apart from the return value
        template<class..., typename Myself = const R&, typename = std:: enable_if_t<has_method_front_ref< Myself >{}>>
        auto operator*() const
            -> decltype(auto) // should return a reference
        {
            static_assert(       has_method_front_ref<Myself>{}  ,"");
            this->m_range_pointer             || [](){throw std:: runtime_error("range:: range-based-for (*end)"    );return false;}();
            (!this->m_range_pointer->empty()) || [](){throw std:: runtime_error("range:: range-based-for (*empty())");return false;}();
            return this->m_range_pointer->front_ref();
        }

        template<class..., typename Myself = const R&, typename = std:: enable_if_t<!has_method_front_ref< Myself >{}>>
        auto operator*() const -> typename R::value_type {
            static_assert(      !has_method_front_ref<Myself>{}  ,"");
            this->m_range_pointer             || [](){throw std:: runtime_error("range:: range-based-for (*end)"    );return false;}();
            (!this->m_range_pointer->empty()) || [](){throw std:: runtime_error("range:: range-based-for (*empty())");return false;}();
            return this->m_range_pointer->front_val();
        }
    };
    template<typename R, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R>{}> >
    auto end(R &) {
        static_assert( !std:: is_reference< R >{}, "");
        static_assert( std::is_base_of<range_tag, R>{}, "" );
        return begin_end_for_range_for<R>{nullptr};
    }
    template<typename R, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R>{}> >
    auto begin(R &r) {
        auto the_end = end(r);
        the_end.m_range_pointer = std:: make_unique<R>( std::move(r) );
        return the_end;
    }

    template<typename R, typename Rnonref, typename F>
    struct filtered_range : private Rnonref {
        static_assert( std:: is_same<Rnonref, std:: remove_reference_t<R>>{} ,"");
        F filter;
        using value_type = typename Rnonref :: value_type;
        static_assert( std:: is_same<value_type, typename Rnonref :: value_type>{} ,""); // just to suppress a gcc warning
        using Rnonref :: empty;
        using Rnonref :: front_val;
        //using Rnonref :: front_ref;
        template<class..., typename U = Rnonref>
        auto front_ref() const -> decltype( U:: front_ref() ) {
            return Rnonref:: front_ref();
        }
        filtered_range(R &&r_, F filter_) : Rnonref( std::forward<R>(r_) ), filter(filter_) {
            skip_unmatching_ones();
        }
        template<class..., typename U = Rnonref>
        auto skip_unmatching_ones() -> std:: enable_if_t< has_method_front_ref< U >{}>
        {
            while(!empty() && !filter(Rnonref:: front_ref())) { // BUG - should make this ref, if possible
                Rnonref:: advance(); // must be done to the base class, as advance will also be replaced
            }
        }
        template<class..., typename U = Rnonref>
        auto skip_unmatching_ones() -> std:: enable_if_t< !has_method_front_ref< U >{}>
        {
            while(!empty() && !filter(Rnonref:: front_val())) { // BUG - should make this ref, if possible
                Rnonref:: advance(); // must be done to the base class, as advance will also be replaced
            }
        }
        void advance() {
            // do at least one advance, but also skipping ones that don't match
            Rnonref :: advance();
            skip_unmatching_ones();
        }
    };
    template<typename R, typename F, class..., typename Rnonref = std::remove_reference_t<R>, typename = std:: enable_if_t<std::is_base_of<range_tag, Rnonref >{}> >
    auto filter(R &&r, F filter) {
        return filtered_range<R, Rnonref, F>{ std::forward<R>(r), filter };
    }

    template<typename R, typename Rnonref, typename F>
    struct transformed_range : private Rnonref {
        static_assert( std:: is_same<Rnonref, std:: remove_reference_t<R>>{} ,"");
        using value_type = decltype(  std:: declval<F>()( std:: declval<typename Rnonref::value_type>() ) );

        F transform;

        transformed_range(R &&r_, F transform_) : Rnonref( std::forward<R>(r_) ), transform(std::move(transform_)) {
        }

        using Rnonref :: empty;
        using Rnonref :: advance;

        value_type front_val() const {
            return transform( Rnonref:: front_val() );
        }
    };
    template<typename R, typename F, class..., typename Rnonref = std::remove_reference_t<R>, typename = std:: enable_if_t<std::is_base_of<range_tag, Rnonref >{}> >
    auto transform(R && r, F && f) {
        return transformed_range<R, Rnonref, F>{ std::forward<R>(r), std:: forward<F>(f) };
    }

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
    template<typename R, class..., typename Rnonref = std::remove_reference_t<R>, typename = std:: enable_if_t<std::is_base_of<range_tag, Rnonref >{}> >
    auto sum (R r) {
        using value_type = typename R :: value_type;
        static_assert( !std::is_same<value_type, bool>{} , "" ); // using bool here is pretty dumb. Usually.
        value_type t = 0;
        for(;!r.empty(); r.advance()) {
            t += r.front_val();
        }
        return t;
    }
    template<typename value_type, typename R, class..., typename Rnonref = std::remove_reference_t<R>, typename = std:: enable_if_t<std::is_base_of<range_tag, Rnonref >{}> >
    auto sum (R r, value_type t) {
        for(;!r.empty(); r.advance()) {
            t += r.front_val();
        }
        return t;
    }

    template<typename R, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R >{}> >
    auto adapt_from_pull(R range_with_pull) {
        struct adaptor_from_pull : public range_tag {
            using value_type = typename R :: value_type;
            R m_range_with_pull;
            value_type m_current_value;
            bool m_is_empty = false;
            adaptor_from_pull(R && range_with_pull) : m_range_with_pull(std:: move(range_with_pull)) {
                advance();
            }
            bool empty() const {
                return m_is_empty;
            }
            void advance() {
                !m_is_empty || [](){throw std:: runtime_error("trying to advance an already-exhausted pull-range");return false;}();
                try {
                    m_current_value = m_range_with_pull.pull();
                } catch (range:: pull_from_empty_range_error &e) {
                    // not sure what start m_current_value will be in.
                    // Don't care either I guess!
                    m_is_empty = true;
                }
            }
            value_type front_val() const {
                !m_is_empty || [](){throw std:: runtime_error("trying to advance an already-exhausted pull-range");return false;}();
                return m_current_value;
            }
            value_type pull() {
                if(empty())
                    throw pull_from_empty_range_error();
                value_type ret = front_val();
                advance();
                return ret;
            }
        };
        return adaptor_from_pull( std:: move(range_with_pull) );
    }


    inline
    auto pull_token_range_from_line(std:: string &full_line, std:: string delims = " \t\n\r") {
        using std:: string;
        using std:: cout;
        using std:: endl;
        struct string_range : public range_tag {
            using value_type = string;
            string const &m_full_line_ref;
            size_t m_pos;
            std:: string m_delims;
            string_range(string &full_line, std:: string delims) :m_full_line_ref(full_line), m_pos(0), m_delims(delims) {
            }
            value_type pull() {
                if(m_pos == string:: npos)
                    throw pull_from_empty_range_error();
                size_t next_white_space = m_full_line_ref.find_first_of(m_delims.c_str(), m_pos);
                string the_token = m_full_line_ref.substr(m_pos, next_white_space-m_pos);
                if(next_white_space == string::npos) {
                    m_pos = string:: npos;
                }
                else {
                    m_pos = next_white_space+1;
                }
                return the_token;
            }
        };
        return string_range(full_line, delims);
    }
    inline
    auto token_range_from_line(std:: string &full_line, std:: string delims = " \t\n\r") {
        return adapt_from_pull(pull_token_range_from_line(full_line, delims));
    }


    template<typename R, class F, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R >{} && has_method_pull<R>{} > >
    auto transform_to_vector(R r, F f)
    ->  std:: vector< decltype(f(r.pull())) >
    {
        std:: vector< decltype(f(r.pull())) > output;
        static_assert( has_method_pull<R>{}, "decltype shoulda caught this");
        for(;;) {
            try {
                output.emplace_back( f(r.pull()) );
            }
            catch (range:: pull_from_empty_range_error &e) {
                return output;
            }
        }
    }
    template<typename R, class F, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R >{} && !has_method_pull<R>{} > >
    auto transform_to_vector(R r, F f) -> std:: vector< decltype(f(r.front_val())) >
    {
        static_assert( !has_method_pull<R>{}, "decltype shoulda caught this");
        std:: vector< decltype(f(r.front_val())) > output;
        while(!r.empty()) {
            output.emplace_back( f(r.front_val()) );
            r.advance();
        }
        return output;
    }
    template<typename R, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R >{} > >
    auto transform_to_vector(R r) {
        return transform_to_vector( std::move(r), [](auto x){return x;} );
    }

    template<typename R, class...
        , typename = std:: enable_if_t<std::is_base_of<range_tag, R >{} >
        , typename = decltype(R :: underlying)
        >
    auto sort(R &r) {
        sort(r.underlying.begin(), r.underlying.end());
    }

    template<typename Source, typename State, typename L>
    auto generic_thing_to__front_range(Source&& t, State && s, L && l) {
        struct gttfr : public range:: range_tag {

            Source m_src;
            State  m_st;
            L      m_l;

            gttfr(Source&&t, State&&s, L&&l) : m_src(std::forward<decltype(t)>(t))
                                             , m_st (std::forward<decltype(s)>(s))
                                             , m_l  (std::forward<decltype(l)>(l)) { }

            using value_type = decltype( m_l(m_src,m_st) );

            static_assert(  !std::  is_reference< value_type > :: value   , "");

            value_type pull()
            { return m_l(m_src,m_st); }
        };
        return gttfr(std::forward<Source>(t), std::forward<State>(s), std::forward<L>(l));
    };
    template<typename T, typename S, typename L>
    auto adapt_pull_function_to_range       (T&& t, S && s, L && l) {
        return adapt_from_pull( generic_thing_to__front_range(std::forward<T>(t), std::forward<S>(s), std::forward<L>(l)) );
    }

    template<typename R, class..., typename = std:: enable_if_t<std::is_base_of<range_tag, R >{} > >
    R copy(R r) {
        return r;
    }

    template<typename T
            , bool is_infinite_ = false
            , class ...
            , typename = std:: enable_if_t<
                                            std:: is_same<T, int64_t>::value
                                          ||std:: is_same<T, int    >::value
                                          ||std:: is_same<T, size_t >::value
                                          ||std:: is_same<T, float>::value
                                          ||std:: is_same<T, double>::value
                                          ||std:: is_same<T, long double>::value
                                          >
            >
    auto range_impl(T b, T e, T step) {
        struct range_impl_struct
            : public range_tag
        {
            using value_type = T;

            struct Members {
                value_type current;
                value_type e;
                value_type step;
            } m_;
            range_impl_struct(Members&& members) : m_( std::move(members) ) {}

            bool empty() const {
                // I guess I should also check for is_infinite here!
                return  m_.current >= m_.e;
            }
            void advance() {
                m_.current += m_.step;
            }
            value_type front_val() const {
                return m_.current;
            }
            static constexpr bool is_infinite() {
                return is_infinite_;
            }
        };
        return range_impl_struct({{b,e,step}});
    }
    template<typename T
            , class ...
            , typename = std:: enable_if_t< std:: is_same<T, int64_t>::value
                                          ||std:: is_same<T, int    >::value
                                          ||std:: is_same<T, size_t >::value
                                          ||std:: is_same<T, float>::value
                                          ||std:: is_same<T, double>::value
                                          ||std:: is_same<T, long double>::value > >
    inline auto iota(T e)
    { return range_impl<T>(0,e,1); }

    template<typename T = int>
    auto infinite_range(T b = 0)
    { return range_impl<T, true>(b, std::numeric_limits<T>::max() ,1); }

    template        <typename R >
    auto max_element(R &&r) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::max_element( std::forward<R>(r).underlying.begin()
                           , std::forward<R>(r).underlying.end()
                           ))
    template        <typename R, typename G >
    auto shuffle(R &&r, G &&g) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::shuffle( std::forward<R>(r).underlying.begin()
                       , std::forward<R>(r).underlying.end()
                       , std::forward<G>(g)
                       ))
    template        <typename R, typename G >
    auto accumulate(R &&r, G &&g) -> AMD_RANGE_DECLTYPE_AND_RETURN(
           std::accumulate( std::forward<R>(r).underlying.begin()
                       , std::forward<R>(r).underlying.end()
                       , std::forward<G>(g)
                       ))
    template        <typename R>
    auto next_permutation(R &&r) -> AMD_RANGE_DECLTYPE_AND_RETURN(
    std::next_permutation( std::forward<R>(r).underlying.begin()
                         , std::forward<R>(r).underlying.end()
                         ))

    // a simpler begin/end thing. Wrote it from scratch within the ssimp project
    template <typename Ib_t, typename Ie_t>
    struct range_from_begin_end_t {
        Ib_t m_b;
        Ie_t m_e;
        auto begin() const { return m_b; }
        auto end  () const { return m_e; }
        bool empty() const { return m_b == m_e; }
        void advance()     { ++m_b; }
        auto current_it() const { return m_b; }
    };

    template <typename Ib_t, typename Ie_t>
    auto range_from_begin_end(Ib_t b, Ie_t e) {
        return range_from_begin_end_t<Ib_t,Ie_t>{move(b),move(e)};
    }
} // namespace range
