/*
 * Aaron McDaid - redoing my range library. Calling it rr.hh for now
 * with namespace 'rr'
 */
#include<utility>
#include<vector>

namespace rr {
    template<typename R, typename = void>
    struct traits;

    namespace impl {
        template<int i>
        struct priority_tag;
        template<int i>
        struct priority_tag : public priority_tag<i-1> {};
        template<>
        struct priority_tag<0> {};
    }

    namespace impl {
        template<typename Possible_Range>
        auto is_range(impl:: priority_tag<2>, Possible_Range &&)
        -> decltype( typename traits< std::remove_reference_t<Possible_Range> > :: value_type{} , std:: true_type{})
        { return {}; }

        template<typename Possible_Range>
        auto is_range(impl:: priority_tag<1>, Possible_Range &&)
        -> std:: false_type
        { return {}; }

        template<typename Possible_Range>
        constexpr
        bool is_range()
        {
            using return_type = decltype(impl:: is_range(impl:: priority_tag<9>{}, std::declval<Possible_Range>()));
            return return_type:: value;
        }
    }
    template<typename Possible_Range>
    constexpr bool is_range_v = impl:: is_range<Possible_Range>();

    template<typename R>
    auto empty  (R const &r)
    ->decltype(traits<R>::empty(r)) {
        return traits<R>::empty(r); }
    template<typename R>
    auto front_val  (R const &r)
    ->decltype(traits<R>::front_val(r)) {
        return traits<R>::front_val(r); }
    template<typename R>
    auto advance    (R       &r)
    ->decltype(traits<R>::advance(r)) {
        return traits<R>::advance(r); }
    template<typename R>
    auto begin      (R       &r)
    ->decltype(traits<R>::begin  (r)) {
        return traits<R>::begin  (r); }
    template<typename R>
    auto end        (R       &r)
    ->decltype(traits<R>::end    (r)) {
        return traits<R>::end    (r); }

    template<typename B, typename E>
    struct pair_of_iterators : public std::pair<B,E>
    {
        static_assert(!std::is_reference<B>{}, "");
        static_assert(!std::is_reference<E>{}, "");
        pair_of_iterators(B b, E e) : std::pair<B,E>(b,e) {}
        /* This struct looks pointless, but it's not.
         * This struct, because it's in the rr:: namespace,
         * can be found by ADL and therefore our begin/end are found easily
         *
         * The main place you see this is in the return from as_range()
         */
    };

    template<typename T>
    struct pair_of_values { T m_begin;
                            T m_end; };

    template<typename I>
    struct iter_is_own_value {
        I m_i;

        bool    operator!=  (iter_is_own_value const & other) const { return  m_i != other.m_i; }
        void    operator++  ()                                      {       ++m_i; }
        I       operator*   ()                                const { return  m_i; }
    };


    template<typename T>
    struct traits< pair_of_values<T> > {
        using R = pair_of_values<T>;
        using value_type = T;
        static
        bool empty      (R const &r) { return r.m_begin == r.m_end ;}
        static
        T    front_val  (R const &r) { return r.m_begin; }
        static
        void advance    (R       &r) {     ++ r.m_begin; }
        static
        auto begin      (R       &r) { return iter_is_own_value<T>{r.m_begin};}
        static
        auto end        (R       &r) { return iter_is_own_value<T>{r.m_end  };}
    };

    inline
    pair_of_values<int> ints(int u) { return {0,u}; }
    inline
    pair_of_values<int> ints(int l, int u) { return {l,u}; }

    template <typename T>
    auto
    as_range(T &v)
    -> pair_of_iterators<   decltype(v.begin())
                        ,   decltype(v.end  ())
                        >
    {
        return {v.begin(),v.end()};
    }

    template <typename T>
    auto
    as_range(T b, T e)
    ->decltype(pair_of_iterators<   decltype(b)
                                ,   decltype(e)
                                >   {b,e})
    {
        return {b,e};
    }

    template<typename I>
    struct traits<std:: pair<I,I>> {
        using R = std:: pair<I,I>;
        using value_type = typename I:: value_type;
        static
        bool empty      (R const &r) {
            return r.first == r.second ;}
        static
        void advance    (R       &r) {
                ++ r.first  ;}
        static
        value_type front_val      (R const &r) {
            return *r.first ;}
    };

    template<typename B, typename E>
    struct traits<pair_of_iterators<B,E>> {
        using R = pair_of_iterators<B,E>;
        using value_type = typename B:: value_type;
        static
        bool empty      (R const &r) {
            return r.first == r.second ;}
        static
        void advance    (R       &r) {
                ++ r.first  ;}
        static
        value_type front_val      (R const &r) {
            return *r.first ;}
        static
        auto begin      (R       &r) { return r.first; }
        static
        auto end        (R       &r) { return r.second; }
    };

    template<typename F, typename Tag_type>
    struct forward_this_with_a_tag {
        F m_r; // may be lvalue or rvalue
        static_assert( is_range_v<F>, "");
    };

    template<typename R, typename F>
    struct mapping_range {
        static_assert(!std::is_reference<R>{},"");
        static_assert(!std::is_reference<F>{},"");
        R m_r;
        F m_f;
    };

    template<typename under_R, typename F>
    struct traits<mapping_range<under_R,F>> {
        using R = mapping_range<under_R,F>;
        using value_type = decltype( rr::front_val  ( std::declval<R>().m_r ));
        static
        bool empty      (R const &r) { return rr:: empty(r.m_r);}
        static
        void advance    (R       &r) { rr::advance( r.m_r ) ;}
        static
        auto front_val      (R const &r) { return r.m_f(rr::front_val  ( r.m_r )) ;}
    };

    template<typename Tag_type>
    struct tagger_t {
    };


    struct map_tag_t            {};     extern  tagger_t<map_tag_t          >   map_range;
    struct map_collect_tag_t    {};     extern  tagger_t<map_collect_tag_t  >   map_collect;
    struct collect_tag_t        {};     extern  collect_tag_t                   collect;    // no need for 'tagger_t', this directly runs
    struct take_collect_tag_t   {};     extern  tagger_t<take_collect_tag_t >   take_collect;

    template<typename R, typename Tag_type
        , std::enable_if_t< is_range_v<R> > * = nullptr
        >
    auto operator| (R && r, tagger_t<Tag_type>) {
        static_assert( is_range_v<R> ,"");
        return forward_this_with_a_tag<R, Tag_type>    {   std::forward<R>(r)  };
    }

    // next, forward 'as_range()' if the lhs is not a range
    template<typename R, typename Tag_type , std::enable_if_t<!is_range_v<R> > * = nullptr >
    auto operator| (R && r, tagger_t<Tag_type> tag) {
        return as_range(std::forward<R>(r)) | tag;
    }

    template<typename R, typename Func>
    auto operator| (forward_this_with_a_tag<R,map_tag_t> f, Func && func) {
        return mapping_range<   std::remove_reference_t<R>
                            ,   std::remove_reference_t<Func>
                            > { std::forward<R   >(f.m_r)
                              , std::forward<Func>(func)
                              };
    }

    template<typename R, typename Func>
    auto operator| (forward_this_with_a_tag<R,map_collect_tag_t> f, Func func) {

        auto r = std::forward<R   >(f.m_r); // copy/move the range here

        using value_type = decltype (   func   (  rr::front_val( r )  ));
        std:: vector<value_type> res;

        while(!rr::empty(r)) {
            res.push_back( func(rr::front_val(r)));
            rr::advance(r);
        }

        return res;
    }

    template<typename R
        , std::enable_if_t< is_range_v<R> > * = nullptr
        >
    auto operator| (R && r, collect_tag_t) {
        static_assert( is_range_v<R> ,"");
        using value_type = decltype (   rr::front_val( r )  );
        std:: vector<value_type> res;

        while(!rr::empty(r)) {
            res.push_back( rr::front_val(r) );
            rr::advance(r);
        }

        return res;
    }

    // next, forward 'as_range()' if the lhs is not a range
    template<typename R , std::enable_if_t< !is_range_v<R> > * = nullptr >
    auto operator| (R && r, collect_tag_t) {
        return as_range(std::forward<R>(r)) | rr:: collect;
    }

    template<typename R>
    auto operator| (forward_this_with_a_tag<R,take_collect_tag_t> f, int how_many) {

        auto r = std::forward<R   >(f.m_r); // copy/move the range here

        using value_type = decltype (   rr::front_val( r )  );
        std:: vector<value_type> res;

        while(!rr::empty(r) && how_many>0) {
            res.push_back( rr::front_val(r));
            rr::advance(r);
            --how_many;
        }

        return res;
    }
#if 0
#endif
} // namespace rr
