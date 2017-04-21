#include"range.hh"
#include<functional>

namespace view {
    struct {} ref_wraps;
    struct {} foreach;
    struct {} unzip_foreach;
    struct {} collect;
    struct {} unzip_collect_transpose;
    struct {} unzip_map;

    template<typename R, typename tag>
    struct temporary_tagged_holder {
        static_assert(!std:: is_reference<R>{}, "");
        R m_r;
        temporary_tagged_holder(R r) : m_r(std::move(r)) {}
    };

    template<typename R>
    auto operator| (R r, decltype(ref_wraps) );

    template<typename R>
    struct ref_wraps_impl {
        R m_range_with_front_ref;

        bool            empty()         const   { return   range:: empty    (m_range_with_front_ref)  ; }
        auto            front_val()     const   { return std:: ref( range:: front_ref(m_range_with_front_ref) ); }
        void            advance()               { range:: advance(m_range_with_front_ref); }
    };

    template<typename R>
    auto operator| (R r, decltype(ref_wraps) ) {
        return ref_wraps_impl<R> {move(r)};
    }

    template<typename R>
    auto operator| (R r, decltype(unzip_foreach) ) {
        return temporary_tagged_holder<R, decltype(unzip_foreach)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(unzip_map) ) {
        return temporary_tagged_holder<R, decltype(unzip_map)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(foreach) ) {
        return temporary_tagged_holder<R, decltype(foreach)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(collect) )
    //-> std:: vector<int>
    {
        using value_type = std:: decay_t< decltype( range:: pull(r) ) >;
        std:: vector<value_type> v;
        while(!r.empty()) {
            v.push_back( range:: pull(r) );
        }
        return v;
    }
    template<size_t ...Is, typename R>
    auto operator_pipe_impl (R r, decltype(unzip_collect_transpose), std:: index_sequence<Is...> )
    //-> std:: vector<int>
    {
        return std::make_tuple( std:: get<Is>( move( r.m_ranges)) | view:: collect ... );
    }
    template<typename R>
    auto operator| (R r, decltype(unzip_collect_transpose) )
    {
        return operator_pipe_impl( move(r), unzip_collect_transpose
                , std:: make_index_sequence< r.width_v >{}
                );
    }

    template<typename R, typename F>
    void operator| (temporary_tagged_holder<R, decltype(unzip_foreach)> r_holder, F && f)
    {
        while(!r_holder.m_r.empty()) {
            utils:: apply   ( std::forward<F>(f)
                            , range:: pull(r_holder.m_r)
            );
        }
    }
    template<typename R, typename F>
    auto operator| (temporary_tagged_holder<R, decltype(unzip_map)> r_holder, F f)
    {
        struct unzip_map_t {
            R m_r;
            F f;
            bool empty() const {
                return m_r.empty();
            }
            auto front_val() {
                return
                utils:: apply   ( std::forward<F>(f)
                                , range:: front_val(m_r)
                                );
            }
            void advance() { m_r.advance(); }
        };
        return unzip_map_t{ move(r_holder.m_r), std:: move(f) };
    }
    template<typename R, typename F>
    void operator| (temporary_tagged_holder<R, decltype(foreach)> r_holder, F && f)
    {
        while(!r_holder.m_r.empty()) {
            std::forward<F>(f)(
                range:: pull(r_holder.m_r)
            );
        }
    }

} // namespace view
