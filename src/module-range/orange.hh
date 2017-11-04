/*
 * orange  - yet another range library
 *
 * By Aaron McDaid - aaron.mcdaid@gmail.com
 *
 *
 * In the code below, I'm trying to organize it so that it can be
 * read from top to bottom and is understandable.
 *
 * Where possible, I put some tests into a 'testing_namespace'. You might
 * find it useful to just search for that string in this file and read
 * the tests.
 *
 *
 * Brief description, and overview of this code
 * ============================================
 *
 * ( This documentation includes some stuff that isn't implemented. We
 *   should implement more! )
 *
 *      vector<int>  v {2,3,5,7};
 *
 *      // print every value
 *      v |foreach| [](auto x) { std::cout << x << '\n'; };
 *
 *      // print the square of each value
 *      v   |mapr|      [](auto x) { return x*x; }
 *          |foreach|   [](auto y) { std:: cout "x^2=" << y << '\n';};
 *
 *      // filter to only the odd ones, then print them:
 *      v   |filter|  [](auto x) { return x%2 == 1; }
 *          |foreach| [](auto x) { std::cout << x << '\n'; };
 *
 * ( say 'mapr' instead of 'map' simply to avoid clashing with 'std::map')
 *
 * Many different types can be considered as 'range types'. The obvious
 * example is a pair of iterators, but there are many others too.
 * A 'vector' is not itself a range; but it is trivially convertable
 * to a range.
 *
 * A 'range' is typically a very lightweight object, like a pointer, that
 * can be copied easily. It doesn't usually "own" the underlying data,
 * but this library supports ownership where appropriate.
 *
 * A range will support some subset of these actions:
 *  -   empty       ::: No more input is available to read
 *  -   front_val   ::: Read the current value - repeated calls will return
 *                      the same value (unless the underlying container has
 *                      been modified by some other part of the system.
 *  -   advance     ::: skip the current item and move to the next
 *
 *  -   full        ::: if an output range can no longer be written to
 *  -   front_ref   ::: return a reference to the current item. Repeated
 *                      calls will return a reference to the same object.
 *  -   push        ::: write a value to an output range and advance. This
 *                      is useful when treating the standard output as an
 *                      output range. It's not possible to define 'front_ref'
 *                      on such a range as we can't meaningful write to the
 *                      "same place" in the output repeatedly. Once we write
 *                      to the stream, our next write must be to the following
 *                      'position' in the output range.
 *  -   pull        ::: return the current value and also advance. As if
 *                      running front_val and then advance. Useful when a
 *                      range doesn't allow repeating read
 *
 * Via traits (see below), you can specify, for your own types, how these
 * actions are to be performed on your objects.
 * These names are all available in the orange:: namespace. They will use
 * the underlying traits actions where they are provided and, in some cases,
 * this library can synethesize extra functions where they are not explicit
 * in your trait; for example, we can synthesize 'orange::pull' from 'front_val'
 * and 'advance' if your trait does not contain 'advance'
 *
 * This becomes useful when you want to combine a range and a function,
 * and create a new range which exposes range where the function has been
 * applied to each element of the underlying range.
 *
 * ==
 *
 * A 'range type', R, here is a type for which the type traits<R> exists.
 * More precisely, traits<R> can also be default constructed. The traits
 * object has no state, its purpose is simply to record the type-specific
 * details. For example, an input range must be able to support the 'orange::empty'
 * function which tells us if more input is available. For a pair of iterators,
 * this means testing the two iterators to see if they are equal. For a file
 * input stream, we test the stream for end-of-file.
 *
 * (I'll try to document the functions in the order they appear below in
 * the code)
 *
 *  is_range_v      ::: constexpr-bool function to test if a given
 *                      type R has a suitable traits<R>
 *
 * The code then has the traits definition for a std::pair of iterators.
 * Traits for other types are specified later in this code, but I brought
 * std::pair to the top as it's simple and helps me to explain this system
 *
 *  template<typename I>
 *  struct traits<std:: pair<I,I>>
 *
 * For now, this just means providing 'empty', 'advance' and 'front_val'
 * In future, some functions for trait a pair as an output range should
 * be added, such as 'full' and 'front_ref'.
 *
 * Next, the functions in 'orange::' are defined, relying on the operations
 * provided in the traits object. For example, this defines 'orange::front_val':
 *
 *  template<typename R>
 *  auto front_val  (R const &r)
 *  ->decltype(traits<R>::front_val(r)) {
 *      return traits<R>::front_val(r); }
 *
 * Another overload of 'orange::front_val' could be provided to synthesize
 * front_val where the traits has 'front_ref', but not 'front_val'.
 *
 */

#include<utility>
#include<vector>

/*
 * orange_utils
 *
 * I define 'is_invokable_v' in this namespace as it's range specific and might be useful elsewhere.
 */
namespace orange_utils {



    /*  priority_tag
     *  ============
     *      'priority_tag' is very useful to specify priority
     *  among overloads that would otherwise be ambiguous.
     *  https://stackoverflow.com/questions/43470741/how-does-eric-nieblers-implementation-of-stdis-function-work
     */
    template<int i>
    struct priority_tag;
    template<int i>
    struct priority_tag : public priority_tag<i-1> {};
    template<>
    struct priority_tag<0> {};


    /*  void_t
     *  ======
     * https://stackoverflow.com/questions/27687389/how-does-void-t-work
     */
    template< typename ... >
    struct voider_t { using type = void; };
    template< typename ... Ts> using void_t = typename voider_t<Ts...> :: type;


    /*  is_invokable_v
     *  ==============
     *  is_invokable_v<F, Args...> tells us if the function object F
     *  can be called with arguments of types Args...
     */
    namespace impl__is_invokable {
        template<typename F, typename ... Args>
        constexpr auto
        is_invokable_one_overload(orange_utils::priority_tag<2>)
        -> decltype( std::declval<F>()(std::declval<Args>()...), true )
        { return true; }

        template<typename F, typename ... Args>
        constexpr auto
        is_invokable_one_overload(orange_utils::priority_tag<1>)
        -> decltype( false )
        { return false; }

        template<typename F, typename ... Args>
        constexpr bool
        is_invokable_v =
                   is_invokable_one_overload<F, Args...>(orange_utils::priority_tag<9>{});
    }

    using impl__is_invokable:: is_invokable_v;  // to 'export' this to the orange_utils namespace


    /* testing_namespace
     * =================
     *  Throughout this file, I'll put tests, using static_assert, into this
     *  namespace. Reading the tests might help you to understand more of
     *  this code.
     */
    namespace testing_namespace {
        /*
         * To make a tester which checks if a give type has a '.size()' method, we define a lambda with the
         * relevant expression 'x.size()'. And also, we test if addition, (x+x), is defined.
         */
        auto checker_for__has_size_method   = [](auto&&x)->decltype(void(  x.size() )){};
        auto checker_for__has_addition      = [](auto&&x)->decltype(void(  x + x    )){};

        template<typename Arg>
        constexpr bool has_size_method  = orange_utils:: is_invokable_v<decltype(checker_for__has_size_method), Arg >;
        template<typename Arg>
        constexpr bool has_addition     = orange_utils:: is_invokable_v<decltype(checker_for__has_addition), Arg >;

        static_assert( has_size_method< std::vector<int> > ,"");
        static_assert(!has_size_method< int              > ,"");
        static_assert( has_size_method< std::vector<int> > ,"");
        static_assert(!has_size_method< int              > ,"");
    }
}

namespace orange {


    /*  traits<R>
     *  =========
     *      If 'R' is a range type, then this traits class tells us
     *  how to use it; how to test if it's empty, for example.
     *  With a pair of iterators, we test for emptiness by testing
     *  if the two iterators equal to each other. With a file input
     *  stream, we would test for emptiness by testing for .eof().
     */
    template<typename R, typename = void> // second template arg is to allow 'void_t' https://stackoverflow.com/questions/27687389/how-does-void-t-work
    struct traits;


    /*  lookup_traits<R>
     *  ================
     *      We don't look up 'traits' directly. We go through 'lookup_traits'
     *  instead, as it drops 'const' and drops references.
     */
    template<typename R
            , typename R_decayed = std::decay_t<R>
            , decltype( traits< R_decayed> {} ) * = nullptr >
    struct lookup_traits : public traits<R_decayed> {};


    /*  checker_for__is_range  is_range_v
     *  =====================  ==========
     *      is_range_v<R> tests if lookup_traits<R> is defined.
     *  This is how we define is a type is a range type or not.
     */
    auto checker_for__is_range=[](auto&&x)->decltype(void(  lookup_traits< decltype(x)>{}  )){};

    template<typename T > constexpr bool
    is_range_v = orange_utils:: is_invokable_v<decltype(checker_for__is_range), T>;


    /*  has_trait_{empty,advance,front_val,front_ref,pull}
     *  ==================================================
     *      In order to 'synthesize' the user-facing functions ( orange::front_val, orange::empty, and so on )
     *  for a range type R, we need a convenient way to check which functions are provided in the trait<R>.
     *  These are the 'has_trait_*' functions defined here:
     */

    auto checker_for__has_trait_empty       = [](auto&&r)->decltype(void( lookup_traits<decltype(r)>::empty    (r) )){};
    auto checker_for__has_trait_advance     = [](auto&&r)->decltype(void( lookup_traits<decltype(r)>::advance  (r) )){};
    auto checker_for__has_trait_front_val   = [](auto&&r)->decltype(void( lookup_traits<decltype(r)>::front_val(r) )){};
    auto checker_for__has_trait_front_ref   = [](auto&&r)->decltype(void( lookup_traits<decltype(r)>::front_ref(r) )){};
    auto checker_for__has_trait_pull        = [](auto&&r)->decltype(void( lookup_traits<decltype(r)>::pull     (r) )){};

    template<typename R> constexpr bool
    has_trait_empty     = orange_utils:: is_invokable_v<decltype(checker_for__has_trait_empty), R>;
    template<typename R> constexpr bool
    has_trait_advance   = orange_utils:: is_invokable_v<decltype(checker_for__has_trait_advance), R>;
    template<typename R> constexpr bool
    has_trait_front_val = orange_utils:: is_invokable_v<decltype(checker_for__has_trait_front_val), R>;
    template<typename R> constexpr bool
    has_trait_front_ref = orange_utils:: is_invokable_v<decltype(checker_for__has_trait_front_ref), R>;
    template<typename R> constexpr bool
    has_trait_pull      = orange_utils:: is_invokable_v<decltype(checker_for__has_trait_pull), R>;


    /*
     * Users will never call the functions in the trait object directly.
     * Instead, we synthesize all the functions, where possible, such
     * as orange:empty, orange::front_val, orange::advance.
     *
     * This design allows us to synthesize extra functions. For example,
     * if a trait has 'front' and 'advance', but not 'pull', then we
     * are still able to synthesize 'orange::pull' using the first two.
     * This allows each trait to focus on the smallest subset of
     * necessary behaviour.
     */


    // just one overload for 'empty'
    template<typename R>
    auto constexpr
    empty  (R const &r)
    ->decltype(lookup_traits<R>::empty(r))
    { return lookup_traits<R>::empty(r); }


    // two overloads for 'front_val', as we can use 'front_ref'
    // instead if it's present.
    template<typename R , std::enable_if_t<
        has_trait_front_val<R&>
    > * = nullptr >
    auto constexpr
    front_val  (R &r)
    ->decltype(auto)
    { return lookup_traits<R>::front_val(r); }

    template<typename R , std::enable_if_t<
        !has_trait_front_val<R&> && has_trait_front_ref<R&>
    > * = nullptr >
    auto constexpr
    front_val  (R &r)
    {   return lookup_traits<R>::front_ref(r); }


    // one overload for 'front_ref'
    template<typename R>
    auto constexpr
    front_ref  (R & r)
    ->decltype(lookup_traits<R>::front_ref(r))
    {   return lookup_traits<R>::front_ref(r); }


    // one overload for 'advance'
    template<typename R>
    auto constexpr
    advance    (R       &r)
    ->decltype(lookup_traits<R>::advance(r))
    {   return lookup_traits<R>::advance(r); }


    /* Next, we see 'begin' and 'end', which are useful
     * for working with range-based for.
     *
     * TODO: synthesize a suitable pair of iterators
     * for range types that don't specify a begin and
     * end of their own.
     */

    // one overload for 'begin'
    template<typename R>
    auto constexpr
    begin      (R       &r)
    ->decltype(lookup_traits<R>::begin  (r))
    {   return lookup_traits<R>::begin  (r); }

    // one overload for 'end'
    template<typename R>
    auto constexpr
    end        (R       &r)
    ->decltype(lookup_traits<R>::end    (r))
    {   return lookup_traits<R>::end    (r); }


    /* Three overloads for 'pull'.
     *  1. has 'pull' in its trait
     *  2. doesn't have 'pull' but does have 'front_val' and 'advance'
     *  3. doesn't have 'pull' nor 'front_val' but does have 'front_ref' and 'advance'
     */
    template<typename R , std::enable_if_t<
        has_trait_pull<R&>
    >* =nullptr>
    auto constexpr
    pull       (R       &r)
    { return lookup_traits<R>::pull     (r); }

    template<typename R , std::enable_if_t<
        !has_trait_pull <R&> && has_trait_front_val<R&> && has_trait_advance<R&>
    >* =nullptr>
    auto constexpr
    pull       (R       &r)
    {
        auto copy = lookup_traits<R>::front_val(r);
        lookup_traits<R>::advance(r);
        return copy;
    }

    template<typename R , std::enable_if_t<
        !has_trait_pull <R&> && !has_trait_front_val<R&> && has_trait_front_ref<R&> && has_trait_advance<R&>
    >* =nullptr>
    auto constexpr
    pull       (R       &r)
    {
        auto copy = lookup_traits<R>::front_ref(r);
        lookup_traits<R>::advance(r);
        return copy;
    }
}



/*
 * Everything above is very general. It's relevant for all range
 * types. Next, we get more specific, by defining some traits
 * and the methods.
 *
 * The functions in the trait are static, and capture the range
 * as R&, where R is a deduced template parameter.
 */

namespace orange {


    // Let's start with the simplest example - a std::pair of iterators
    template<typename I, typename J>
    struct traits<std:: pair<I,J>> {

        template<typename R> static constexpr
        bool
        empty           (R & r)   { return r.first == r.second ;}

        template<typename R> static
        void
        advance         (R & r)   { ++ r.first  ;}

        template<typename R> static constexpr
        decltype(auto)
        front_ref       (R & r)   { return * r.first ;}
    };


    namespace testing_namespace {
        static_assert(is_range_v< std::pair< std::vector<int>::iterator,  std::vector<int>::iterator> >, "");
        static_assert(is_range_v< std::pair<int*, int*> >, "");
        static_assert( has_trait_empty    < std::pair<int*, int*> > , "");
        static_assert(!has_trait_front_val< std::pair<int*, int*> > , "");
        static_assert( has_trait_front_ref< std::pair<int*, int*> > , "");
        static_assert(!has_trait_front_ref< std::vector<int> > , "");
    }


    struct orange_use_the_methods{}; /* if a class inherits from this, then it means
                                      * that it has methods instead of having to
                                      * manually specify the traits */
    template<typename T>
    struct traits< T , orange_utils:: void_t< std::enable_if_t< std::is_base_of<orange_use_the_methods, T>{} > > > {
        static_assert(!std::is_const<T>{} ,"");
        template<typename R> static constexpr
        auto
        empty      (R &  r)
        ->decltype(r.empty())
        {   return r.empty(); }

        template<typename R> static constexpr
        auto
        front_ref  (R &  r)
        ->decltype(r.front_ref())
        {   return r.front_ref(); }

        template<typename R> static constexpr
        auto
        front_val  (R &  r)
        ->decltype(r.front_val())
        {   return r.front_val(); }

        template<typename R> static constexpr
        auto
        pull       (R &  r)
        ->decltype(r.pull     ())
        {   return r.pull     (); }

        template<typename R> static constexpr
        auto
        advance    (R &  r)
        ->decltype(r.advance  ())
        {   return r.advance  (); }
    };

    struct orange_traits_are_static_here{};
    template<typename T>
    struct traits< T , orange_utils:: void_t< std:: enable_if_t<
        std:: is_same<
            typename T::orange_traits_are_static_here
            , orange  ::orange_traits_are_static_here
        >{}
    > >> {


        template<typename R> static constexpr auto
        empty      (R &  r)
        ->decltype(R:: orange_empty    (r))
        {   return R:: orange_empty    (r); }

        template<typename R> static constexpr auto
        advance    (R &  r)
        ->decltype(R:: orange_advance  (r))
        {   return R:: orange_advance  (r); }

        template<typename R> static constexpr auto
        front_ref  (R &  r)
        ->decltype(R:: orange_front_ref(r))
        {   return R:: orange_front_ref(r); }

        template<typename R> static constexpr auto
        front_val  (R &  r)
        ->decltype(R:: orange_front_val(r))
        {   return R:: orange_front_val(r); }

        template<typename R> static constexpr auto
        pull       (R &  r)
        ->decltype(R:: orange_pull     (r))
        {   return R:: orange_pull     (r); }
    };
}

namespace orange {

    /*
     * Next, a 'pair_of_iterators' type in the orange:: namespace. The main (only?)
     * reason for this (as opposed to an std::pair of iterators) is to allow
     * 'begin' and 'end' to be defined appropriately, allowing  for(auto x : r) to
     * work.  This is the class used when applying thing like '|'
     */

    template<typename B, typename E>
    struct pair_of_iterators : public std::pair<B,E>
    {
        static_assert(!std::is_reference<B>{}, "");
        static_assert(!std::is_reference<E>{}, "");

        using std:: pair<B,E> :: pair; // to inherit the constructors

        /* This struct looks pointless, but it's not.
         * This struct, because it's in the orange:: namespace,
         * can be found by ADL and therefore our begin/end are found easily
         *
         * The main place you see this is in the return from as_range()
         */
    };

    /*
     * 'pair_of_values', so that we can range between a pair of numbers.  See
     * the 'ints' function below.
     *
     * This is related to 'iter_is_own_value', which is what we get if we call 'begin'
     * and 'end' on a 'pair_of_values'.
     */
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
        template<typename R> static constexpr
        bool empty      (R &  r)   { return r.m_begin == r.m_end ;}

        template<typename R> static constexpr
        T    front_val  (R &  r)   { return r.m_begin; }

        template<typename R> static constexpr
        void advance    (R &  r)   {     ++ r.m_begin; }

        template<typename R> static constexpr
        auto begin      (R &  r)   { return iter_is_own_value<T>{r.m_begin};}

        template<typename R> static constexpr
        auto end        (R &  r)   { return iter_is_own_value<T>{r.m_end  };}
    };

    inline
    constexpr
    pair_of_values<int> ints(int u) { return {0,u}; }
    inline
    constexpr
    pair_of_values<int> ints(int l, int u) { return {l,u}; }

    /*
     * as_range
     * Converts a non-range to a range, where appropriate. The obvious examples
     * are a container such as 'std::vector' or 'std::list'.  as_range can also
     * be called with two iterators.
     * I should also define an 'as_range' overload that accepts a range and
     * forwards it as-is.
     */
    template <typename T>
    auto constexpr
    as_range(T &v)
    -> pair_of_iterators<   decltype(v.begin()) ,   decltype(v.end  ()) >
    { return {v.begin(),v.end()}; }

    template <typename T, std:: size_t N>
    auto constexpr
    as_range(T (&v)[N])
    -> pair_of_iterators< T*, T* >
    { return {std::begin(v),std::end(v)}; }

    template <typename T>
    auto constexpr
    as_range(T b, T e)
    ->decltype(pair_of_iterators<   decltype(b) ,   decltype(e) >   {b,e})
    { return {b,e}; }

    template<typename B, typename E>
    struct traits<pair_of_iterators<B,E>> {
        template<typename R> static constexpr
        bool
        empty           (R & r)   { return r.first == r.second ;}

        template<typename R> static constexpr
        void
        advance         (R & r)   { ++ r.first  ;}

        template<typename R> static constexpr
        auto
        front_val       (R & r)   { return * std::forward<R>(r) .first ;}

        template<typename R> static constexpr
        decltype(auto)
        front_ref       (R & r)   { return * std::forward<R>(r) .first ;}

        template<typename R> static constexpr
        auto begin      (R & r)   { return r.first; }

        template<typename R> static constexpr
        auto end        (R & r)   { return r.second; }
    };

    template<typename C>
    struct owning_range { // non-copyable
        static_assert(!std::is_reference<C>{}   ,"");
        static_assert(!is_range_v<C>            ,"");

        using R = decltype( as_range(std::declval<C&>()) );
        static_assert(!std::is_reference<R>{} ,"");
        static_assert( is_range_v<R>            ,"");

        C m_c;
        R m_r;


        constexpr
        owning_range    (C && c) // takes only an rvalue reference, *not* a universal reference
        : m_c(std::move(c)) , m_r( as_range(m_c) )
        {}

        // don't allow this to be copied
        owning_range    (owning_range const &) = delete;
        owning_range &  operator=  (owning_range const &) = delete;

        // ... but allow moving
        constexpr
        owning_range    (owning_range      &&) = default;
        constexpr
        owning_range &  operator=  (owning_range      &&) = default;
    };

    template<typename T>
    struct traits< owning_range<T> > {
        template<typename R> static constexpr
        bool empty      (R &  r)   { return orange::empty(r.m_r) ;}
        template<typename R> static constexpr
        auto pull       (R &  r)   { return orange::pull (r.m_r) ;}
    };
    template <typename T
        , std::enable_if_t< !std::is_reference<T>{} > * = nullptr
        >
    auto constexpr
    as_range(T &&t) // enable_if is used here, to ensure it is only an rvalue
    {
        return owning_range<T>{ std::forward<T>(t) };
    }

    /*
     * Above, all the basic underlying technology for a range has
     * been defined. Now, the 'user-facing' code must be implemented,
     * allowing   |mapr|  and  |filter|  and so on.
     *
     *  v |foreach| [](auto x){ std::cout << x << '\n'; }
     *
     * The above works because we overload the '|' operator. 'foreach'
     * is an object of an empty tag type. Therefore  v|foreach  is
     * a valid expression which doesn't do much except capture
     * a copy of the range object
     * Then, via another overload of  |  , we apply the lambda.
     * So, the above can be read as
     *
     *  (v | foreach)  |  [](auto x){ std::cout << x << '\n'; }
     */

    template<typename Tag_type>
    struct tagger_t {
        constexpr tagger_t() {} // clang-3.8.0 insists on a user-provided default constructor
    };

    struct foreach_tag_t        {};     constexpr   tagger_t<foreach_tag_t      >   foreach;
    struct filter_tag_t         {};     constexpr   tagger_t<filter_tag_t       >   filter;
    struct map_collect_tag_t    {};     constexpr   tagger_t<map_collect_tag_t  >   map_collect;
    struct take_collect_tag_t   {};     constexpr   tagger_t<take_collect_tag_t >   take_collect;
    struct map_tag_t            {};     constexpr   tagger_t<map_tag_t          >   map_range;
                                        constexpr   tagger_t<map_tag_t          >   mapr;
    struct collect_tag_t{constexpr collect_tag_t(){}};
                                        constexpr            collect_tag_t          collect;    // no need for 'tagger_t', this directly runs
    struct accumulate_tag_t{constexpr accumulate_tag_t(){}};
                                        constexpr            accumulate_tag_t       accumulate;    // no need for 'tagger_t', this directly runs


    // the type to capture the value, i.e. for the left-hand '|'
    // of   (x|operation|func)
    template<typename R, typename Tag_type>
    struct forward_this_with_a_tag {
        R m_r;
        static_assert(!std:: is_reference<R>{}, "");
        static_assert( is_range_v< R >, "");
    };

    /*
     * The actual overloads of '|' are here.
     */
    template<typename R, typename Tag_type
        , std::enable_if_t< is_range_v<R> > * = nullptr // if 'r' is a range
        >
    auto constexpr
    operator| (R r, tagger_t<Tag_type>) {
        static_assert( is_range_v<R> ,"");
        return forward_this_with_a_tag<R, Tag_type>    {   std::move(r)  };
    }
    template<typename R, typename Tag_type
        , typename Rnonref = std::remove_reference_t<R>
        , std::enable_if_t<!is_range_v<Rnonref> > * = nullptr // if 'nr' is a not a range
        >
    auto constexpr
    operator| (R && nr, tagger_t<Tag_type> tag)
    ->decltype(as_range(std::forward<R>(nr)) | tag)
    {
        return as_range(std::forward<R>(nr)) | tag;
    }

    /*
     * Now, to start defining the various  |operations|
     */

    template<typename R, typename F>
    struct mapping_range {
        static_assert(!std::is_reference<R>{},"");
        static_assert(!std::is_reference<F>{},"");
        static_assert( is_range_v<R>, "");
        R m_r;
        F m_f;

        using orange_traits_are_static_here = orange:: orange_traits_are_static_here;
        template<typename M> static constexpr bool
        orange_empty      (M &m) { return orange:: empty(m.m_r);}
        template<typename M> static constexpr void
        orange_advance    (M &m) { orange::advance( m.m_r ) ;}
        template<typename M> static constexpr auto
        orange_front_val  (M &m) { return m.m_f(orange::front_val  ( m.m_r )) ;}
    };

    // |mapr| or |map_range|
    template<typename R, typename Func>
    auto constexpr
    operator| (forward_this_with_a_tag<R,map_tag_t> f, Func && func) {
        return mapping_range<   std::remove_reference_t<R>      // so we store it by value
                            ,   std::remove_reference_t<Func>
                            > { std::move         (f.m_r)
                              , std::forward<Func>(func)
                              };
    }

    template<typename R, typename F>
    struct filter_range
    : public orange_use_the_methods
    {
        static_assert(!std::is_reference<R>{},"");
        static_assert(!std::is_reference<F>{},"");
        static_assert( is_range_v<R>, "");

        R m_r;
        F m_f;

        constexpr
        void
        skip_if_necessary() {
            while(!orange::empty(m_r) && !m_f(orange::front_val(m_r)))
            { orange::advance(m_r); }
        }

        template<typename RR, typename FF>
        constexpr
        filter_range(RR && r, FF && f)
        : m_r(std::forward<RR>(r)) , m_f(std::forward<FF>(f))
        { skip_if_necessary(); }

        constexpr bool
        empty       () const    { return orange:: empty(m_r); }
        constexpr void
        advance     ()          { orange::advance( m_r ); this->skip_if_necessary() ;}
        constexpr auto
        front_val   () const    { return orange::front_val  ( m_r ) ;}
    };

    // |filter|
    template<typename R, typename Func>
    auto constexpr
    operator| (forward_this_with_a_tag<R,filter_tag_t> f, Func && func) {
        return filter_range <   std::remove_reference_t<R>      // so we store it by value
                            ,   std::remove_reference_t<Func>
                            > { std::move         (f.m_r)
                              , std::forward<Func>(func)
                              };
    }

    // |collect|
    template<typename R, typename Func>
    auto constexpr
    operator| (forward_this_with_a_tag<R,map_collect_tag_t> f, Func func) {

        static_assert(!std::is_reference<decltype(f.m_r)>{}, "");

        using value_type = decltype (   func   (  orange::front_val( f.m_r )  ));
        static_assert(!std::is_reference<value_type>{} ,"");
        std:: vector<value_type> res;

        while(!orange::empty(f.m_r)) {
            res.push_back( func(orange::pull(f.m_r)));
        }

        return res;
    }

    template<typename R
        , typename Rnonref = std::remove_reference_t<R>
        , std::enable_if_t< is_range_v<Rnonref> > * = nullptr
        >
    auto constexpr
    operator| (R r, collect_tag_t) {
        static_assert( is_range_v<R> ,"");
        using value_type = decltype (   orange::front_val( r )  );
        static_assert(!std::is_reference<value_type>{} ,"");
        std:: vector<value_type> res;

        while(!orange::empty(r)) {
            res.push_back( orange::pull(r) );
        }

        return res;
    }

    // next, forward 'as_range()' if the lhs is not a range
    template<typename R
        , typename Rnonref = std::remove_reference_t<R>
        , std::enable_if_t< !is_range_v<Rnonref> > * = nullptr >
    auto constexpr
    operator| (R && r, collect_tag_t) {
        return as_range(std::forward<R>(r)) | orange:: collect;
    }

    //  |accumulate
    template<typename R
        , std::enable_if_t< is_range_v<R> > * = nullptr
        >
    auto constexpr
    operator| (R r, accumulate_tag_t) {
        static_assert(!std::is_reference<R>{},"");
        static_assert( is_range_v<R> ,"");

        using value_type = std::remove_reference_t<decltype(orange::pull(r))>;
        value_type total = 0;

        while(!orange::empty(r)) {
            total += orange::pull(r);
        }

        return total;
    }

    // |take_collect
    template<typename R>
    auto constexpr
    operator| (forward_this_with_a_tag<R,take_collect_tag_t> f, int how_many) {

        static_assert(!std::is_reference<decltype(f.m_r)>{}, "");

        using value_type = decltype (   orange::front_val( f.m_r )  );
        std:: vector<value_type> res;

        while(!orange::empty(f.m_r) && how_many>0) {
            res.push_back( orange::front_val(f.m_r));
            orange::advance(f.m_r);
            --how_many;
        }

        return res;
    }

    namespace testing_namespace {
        static_assert( 10 ==  (ints(5) | accumulate)  ,"");
        constexpr double x[] = {1.0, 2.7, 3.14};
        static_assert(1.0 + 2.7 + 3.14 == (as_range(std::begin(x), std::end(x)) | accumulate) ,"");
        static_assert(1.0 + 2.7 + 3.14 == (as_range(x)                          | accumulate) ,"");

        struct greater_than_5_t {
            constexpr greater_than_5_t() {}
            constexpr bool operator() (int x) const { return x>5; }
        };
        struct odd_t {
            constexpr odd_t() {}
            constexpr bool operator() (int x) const { return x % 2 == 1; }
        };
        struct even_t {
            constexpr even_t() {}
            constexpr bool operator() (int x) const { return x % 2 == 0; }
        };
        struct negate_t {
            constexpr negate_t() {}
            template<typename T>
            constexpr auto operator() (T x) const { return -x; }
        };
        static_assert(20 == (ints(10) |filter| even_t{}             |accumulate) ,"");
        static_assert(25 == (ints(10) |filter| odd_t{}              |accumulate) ,"");
        static_assert(30 == (ints(10) |filter| greater_than_5_t{}   |accumulate) ,"");
        static_assert(-30 == (ints(10) |filter| greater_than_5_t{} |mapr| negate_t{}   |accumulate) ,"");

        template<typename T, std::size_t N>
        struct orange_over_an_array
        {   // This class is only useful in the constexpr tests. It's rubbish, but it plays nicely with
            // constexpr. It's the only way I can think of to get generic data into a range.
            T m_array[N];
            using Tc = T const;
            std:: size_t offset;

            constexpr
            auto empty() const { return offset >= N; }
            constexpr
            auto pull ()       { return m_array[offset++]; }
            constexpr
            T & front_ref ()       { return m_array[offset]; }
            constexpr
            T   front_val ()       { return m_array[offset]; }
        };
    }
    // (dropping back up a namespace temporarily, just
    // to define the trait for 'orange_over_an_array'
    template<typename T, std::size_t N>
    struct traits<testing_namespace:: orange_over_an_array<T,N>> {

        template<typename RR>
        static constexpr
        decltype(auto) empty(RR & r) { return r.empty(); }

        template<typename RR>
        static constexpr
        decltype(auto) pull(RR & r) { return r.pull(); }

        template<typename RR>
        static constexpr
        decltype(auto)
        front_ref(RR & r)
        {   return r.front_ref(); }
    };
    namespace testing_namespace{
        static_assert(60 == (orange:: testing_namespace:: orange_over_an_array<int, 3>({{10,20,30},0}) | accumulate) ,"");

        constexpr
        int foo() {
            auto ooaa = orange_over_an_array<int, 3>({{10,20,30},0});
            orange:: front_ref(ooaa) += 100;
            return ooaa | accumulate;
        }
        static_assert( 160 == foo() , "");
    }
} // namespace orange
