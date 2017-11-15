/* A simple string formatting library, loosely based on that of C#, but using
 * a compile-time string to be type-safe and efficient.
 *
 * I'll start with a simple run-time implementation first, as it's easier
 * to develop and test.
 *
 * Based on C#: https://msdn.microsoft.com/en-us/library/txafckwd(v=vs.110).aspx
 *
 * { index,alignment;delimiter:formatString}
 *
 * alignment:
 *      - the number is a minimum width, e.g. '>5'
 *      - ?consider '>' to specify a maximum width?
 *      - ?consider '=' to specify an exact width?
 *      - negative means left aligned
 *
 * delimiter
 */
#ifndef HH_FORMAT_HH
#define HH_FORMAT_HH

#include<string>
#include<cstring>
#include<cassert>
#include<vector>
#include<memory>

#include "../bits.and.pieces/utils.hh"
namespace utils {
    template<typename T, T c>
    constexpr compile_time_constant_as_a_type<T, c>     cx_val  = compile_time_constant_as_a_type<T,c>{}; // Consider using my own type here instead of integral_constant?
} // namespace utils

#define AMD_FORMATTED_STRING(s, ...) format:: do_formatting( [](){ struct local { constexpr static char const * str() { return s; }}; return format:: make_a_char_pack_from_stringy_type<local>::type{}; }() ,__VA_ARGS__)

namespace format {
    struct string_view {
        char const *    backing;
        size_t          b;
        size_t          e;
        explicit
        string_view     (char const *s)
            :   backing(s)
            ,   b(0)
            ,   e(std:: strlen(s))
        { }
        std:: string            get_as_string()     const {
            assert(b <= e);
            assert(e <= std::strlen(backing));
            return std::string(backing + b, e-b);
        }
        char                    pop_front() {
            assert(!empty());
            char c = at(0);
            ++b;
            return c;
        }
        bool                    empty()     const {
            assert(b<=e);
            return b==e;
        }
        char                    at(size_t i) const {
            assert(b+i < e);
            return backing[b+i];
        }
        char                    sneak_ahead() {
            assert(backing[e] != '\0');
            char snuck_in = backing[e];
            ++e;
            return snuck_in;
        }
        size_t                  size() const {
            return e-b;
        }

        template<size_t Nplus1>
        bool                    operator==(char const (&a)[Nplus1]) {
            assert(a[Nplus1-1] == '\0');
            if(size() != Nplus1-1)
                return false;
            return
                utils:: make_a_pack_and_apply_it<Nplus1-1, size_t> (
                    [&](auto ... idxs) {
                        bool all_tested = utils:: and_all( a[idxs] == backing[b+idxs] ...);
                        return all_tested;
                    }
                );
        }
    };

    inline
    std:: pair<string_view, string_view>  consume_simple(string_view all)     {
        /* Consume one of four things:
         *  - literal '{{', representing a single '{'
         *  - literal '}}', representing a single '}'
         *  - a format string, beginning with '{'
         *  - as many non-{} characters as possible
         */
        string_view     left (all);
        string_view     right(all);
        left.e = left.b;
        assert  (   all.get_as_string()
                ==  left.get_as_string()
                +   right.get_as_string()
                );

        if  (   right.size() >= 2
             && right.at(0) == right.at(1)
             && (   right.at(0) == '{' || right.at(0) == '}'    )
            ) {
            right.pop_front(); // skip over one of them
            char popped = right.pop_front();
            char snuck_in = left.sneak_ahead();
            assert(popped == snuck_in); // == '{' or '}'
            return {left,right};
        }

        if  (   right.size() >=2
             && right.at(0) == '{'
             && right.at(1) != '{'
            ) {
            // a formatting substring - must read up until
            // the corresponding '}'
            int current_stack = 0;
            do {
                if(right.at(0) == '{') {
                    ++current_stack;
                }
                if(right.at(0) == '}') {
                    --current_stack;
                }
                char popped = right.pop_front();
                char snuck_in = left.sneak_ahead();
                assert(popped == snuck_in);
            } while(current_stack > 0);
            assert(current_stack == 0);
            return {left,right};
        }

        // move everything from the beginning of `right` to the end of
        // `left` if it's not a special character
        while(1) {
            if  (   !right.empty()
                 && right.at(0) != '{'
                 && right.at(0) != '}'
                 ) {
                char popped = right.pop_front();
                char snuck_in = left.sneak_ahead();
                assert(popped == snuck_in);
            }
            else {
                break;
            }
        }

        assert  (   all.get_as_string()
                ==  left.get_as_string()
                +   right.get_as_string()
                );
        return {left, right};
    }
    struct stored_format_data_for_later_I {
        virtual     std:: string apply_this_format(string_view) = 0;
    };
    template<typename T>
    struct stored_format_data_for_later : stored_format_data_for_later_I {
        T   m_data;

        template<typename U>
        stored_format_data_for_later(U &&u) : m_data( std::forward<U>(u) ) {}

        virtual     std:: string apply_this_format(string_view) override {
            return "?";
        }
    };

    template<typename ...Ts>
    std:: string    format(char const *fmt, Ts const & ...ts) {
        std:: vector< std:: shared_ptr< stored_format_data_for_later_I >> all_arg_data{
            std:: make_shared<stored_format_data_for_later<Ts>>(ts) ... };

        string_view remainder(fmt);

        while(!remainder.empty()) {
            auto p = consume_simple(remainder);

            remainder = p.second;

            auto current_token = p.first;
            if(current_token == "{0}") {
            }
            if(current_token == "{1}") {
            }
        }
        return fmt;
    }

    /* Start the 'compile-time' version of this */
    constexpr static size_t cx_strlen(char const *s) {
        int l = 0;
        while(*s++ != '\0') {
            ++l;
        }
        return l;
        //return *s=='\0' ? 0 : (1+cx_strlen(s+1));
    }

    template<typename string_provider_t>
    struct make_a_char_pack_from_stringy_type {
        constexpr static size_t len =  cx_strlen( string_provider_t :: str());
        constexpr static char   at(size_t i) {
            return string_provider_t:: str() [i];
        }
        static
        auto compute_the_char_pack_type () {
            return
            utils:: make_a_pack_and_apply_it<len, size_t> ( [](auto ...idxs) {  using return_type = utils:: char_pack< make_a_char_pack_from_stringy_type::at(idxs) ...>; return return_type{}; });
        }
        using type = decltype( compute_the_char_pack_type() );
    };

#define FORMAT_ENABLE_IF_THINGY(...) std:: enable_if_t< __VA_ARGS__ > * = nullptr

    namespace type_vector_ns {
        template<typename ...T>
        struct type_vector {
            using as_tuple = std:: tuple<T...>;

            template<typename V>
            constexpr static auto at_ty(V v) {
                return typename std:: tuple_element< v , as_tuple> :: type {};
                //return std::get<V>(as_tuple{});
            }

            constexpr static size_t     size()  { return sizeof...(T); }
            constexpr static auto       size_ty()   { return utils:: cx_val<size_t,sizeof...(T)>; }
        };

        // make_type_vector
        template<typename ...T>
        type_vector<T...> make_type_vector( T ... ) {
            static_assert( utils:: and_all( true, std:: is_empty<T> ::value ... ) ,"");
            return {};
        }

        // map_type
        template<typename ...T, typename F>
        auto map_type   (   type_vector<T...>   , F f)
        -> decltype(    make_type_vector( f( T{} ) ... )    )
        { return {}; (void)f; }

        // concat_type_vectors
        template<typename ...T1, typename ...T2>
        auto concat_type_vectors   (   type_vector<T1...>   , type_vector<T2...> ) {
            return make_type_vector( T1{}..., T2{}... );
        }

        // cumsum_type
        template<typename I>
        auto cumsum_type   (    I     ,   type_vector<> emp) {
            return emp;
        }
        template<typename T0, typename ...Trest, typename I>
        auto cumsum_type   (    I init,   type_vector<T0, Trest...> tv) {
            auto first_sum = init + T0{};
            auto ret= concat_type_vectors   (   make_type_vector(   first_sum )
                                            ,   cumsum_type     (   first_sum, type_vector<Trest...>{})     );
            static_assert( tv.size() == ret.size() ,"");
            return ret;
        }

        // which
        namespace detail {
            template<size_t off>
            auto    which(type_vector<>) {
                return make_type_vector();
            }
            template<size_t off, typename ...Brest> auto    which(type_vector< utils:: compile_time_constant_as_a_type<bool, true >, Brest...>);
            template<size_t off, typename ...Brest> auto    which(type_vector< utils:: compile_time_constant_as_a_type<bool, false>, Brest...>);
            template<size_t off, typename ...Brest>
            auto    which(type_vector< utils:: compile_time_constant_as_a_type<bool, true>, Brest...>) {
                return concat_type_vectors  (   make_type_vector( utils:: cx_val<size_t, off> )
                                            ,   detail:: which<off+1>( type_vector<Brest...>{} )        );
            }
            template<size_t off, typename ...Brest>
            auto    which(type_vector< utils:: compile_time_constant_as_a_type<bool, false>, Brest...>) {
                return detail:: which<off+1>( type_vector<Brest...>{} );
            }
        }
        template<typename ...B>
        auto    which(type_vector<B...> bvt) {
            return detail:: which<0>(bvt);
        }
    } // namespace type_vector_ns


    template<typename corresponding_char_pack>
    struct plain_output { };

    template<typename corresponding_char_pack>
    struct formatter {
        template<typename Tup
            ,   class ...
            ,   typename copy_of_corresponding_char_pack = corresponding_char_pack
            ,   FORMAT_ENABLE_IF_THINGY( copy_of_corresponding_char_pack::at(0) == '{' )
            ,   FORMAT_ENABLE_IF_THINGY( copy_of_corresponding_char_pack::at(1) >= '0' )
            ,   FORMAT_ENABLE_IF_THINGY( copy_of_corresponding_char_pack::at(1) <= '9' )
            ,   FORMAT_ENABLE_IF_THINGY( copy_of_corresponding_char_pack::at(2) == '}' )    >
        void run(std:: ostringstream &oss, Tup && tup) {
            oss << std::get< copy_of_corresponding_char_pack::at(1) - '0' > (tup);
        }
    };

    template < char ...c >
    auto parse_one_thing(utils:: char_pack<'{', '{', c...> s) {
        auto    head    = s.template substr<1>();
        auto    tail    = s.template substr<2, s.size()>();
        return std:: make_pair( plain_output<decltype(head)>{}, tail);
    }
    template < char ...c >
    auto parse_one_thing(utils:: char_pack<'}', '}', c...> s) {
        auto    head    = s.template substr<1>();
        auto    tail    = s.template substr<2, s.size()>();
        return std:: make_pair( plain_output<decltype(head)>{}, tail);
    }

    template <char first_char, char ...c
        , typename ...
        , FORMAT_ENABLE_IF_THINGY( first_char != '{' && first_char != '}')
        >
    auto parse_one_thing(utils:: char_pack<first_char, c...> s) {
        // parse everything up to, but not including, '{' or '}' or '\0'
        size_t constexpr length_of_head = s.find_first_of('{','}','\0');
        static_assert( length_of_head > 0 ,"");

        auto    head    = s.template substr<length_of_head>();
        auto    tail    = s.template substr<length_of_head, s.size()>();

        return std:: make_pair( plain_output<decltype(head)>{}, tail);
    }

    template <char m, char ...c
        , FORMAT_ENABLE_IF_THINGY( m != '{' )
        >
    auto parse_one_thing(utils:: char_pack<'{', m, '}', c...> s) {
        using utils:: make_a_pack_and_apply_it;
        using utils:: cx_val;
        using type_vector_ns:: make_type_vector;
        using type_vector_ns:: which;

        auto chars_as_type_vector = make_a_pack_and_apply_it<s.size(), size_t>(
                [&](auto ... idxs) { return make_type_vector( cx_val<char, s.at(idxs)> ... ) ; });
        auto mapped = map_type(chars_as_type_vector, [](auto x){
                return cx_val<int,      x=='{'  ?    1  :
                                        x=='}'  ?   -1  :
                                                     0       >;});
        auto cumulative_sum = cumsum_type(cx_val<int,0>, mapped);
        auto is_zero        = map_type  (   cumulative_sum
                                        ,   [](auto x) { return x == cx_val<int, 0>; });
        auto which_are_zero = which(is_zero);
        static_assert( which_are_zero.size_ty() > 0, "");
        auto constexpr first_zero = which_are_zero . at_ty( cx_val<int, 0> );

        // Next few asserts are too conservative
        static_assert( first_zero == cx_val<size_t,2> ,""); // first time the braces are balanced

        static_assert( cumulative_sum. at_ty(cx_val<int, 0>) == 1 ,"");
        static_assert( cumulative_sum. at_ty(cx_val<int, 1>) == 1 ,"");
        static_assert( cumulative_sum. at_ty(cx_val<int, 2>) == 0 ,""); // it's the closing brace

        size_t constexpr length_of_head = first_zero+1; // to include the closing brace

        auto    head    = s.template substr<length_of_head>();
        auto    tail    = s.template substr<length_of_head, s.size()>();

        return std:: make_pair( formatter<decltype(head)>{}, tail);
    }

    inline
    auto parse_many_things(utils:: char_pack<>) {
        return std:: make_tuple();
    }
    template <char ...c>
    auto parse_many_things(utils:: char_pack<c...> s) {
        static_assert(sizeof...(c) > 0 ,"");
        auto head_and_tail = parse_one_thing(s);
        return std:: make_pair  (                        head_and_tail.first
                                ,   parse_many_things   (head_and_tail.second)
                );
    }

    template<   size_t          num_formatters_so_far
            ,   typename        Tup         >
    void go_forth_and_print(std::ostringstream &, std::tuple<> , Tup &&) {
    }

    template<   size_t          num_formatters_so_far
            ,   typename        head_subtype
            ,   typename        tail_type
            ,   typename        Tup         >
    void go_forth_and_print(std::ostringstream &oss, std:: pair<plain_output<head_subtype>, tail_type> p, Tup && tup ) {
        oss << head_subtype :: c_str0();
        go_forth_and_print<num_formatters_so_far>(oss, p.second, std::forward<Tup>(tup));
    }

    template<   size_t          num_formatters_so_far
            ,   typename        head_subtype
            ,   typename        tail_type
            ,   typename        Tup         >
    void go_forth_and_print(std::ostringstream &oss, std:: pair<formatter<head_subtype>, tail_type> p, Tup && tup  ) {
        p.first.run(oss, std::forward<Tup>(tup) );
        go_forth_and_print<1+num_formatters_so_far>(oss, p.second, std::forward<Tup>(tup));
    }

    template<char ...chars, typename ...Ts>
    std:: string do_formatting( utils:: char_pack<chars...> s, Ts && ... ts) {

        auto all_parsed = parse_many_things(s);

        std:: ostringstream oss;
        go_forth_and_print<0>(oss, all_parsed, std:: forward_as_tuple(std::forward<Ts>(ts)...) );

        return oss.str();
    }

} // namespace format
#endif
