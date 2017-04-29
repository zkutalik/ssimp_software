#ifndef AMD_RANGE_VIEW_HH
#define AMD_RANGE_VIEW_HH

#include"range.hh"
#include<functional>

namespace range {
namespace view {
    template<typename V>
    auto enumerate_vector(V & v) {
        return zip_val  (range:: ints(v.size())
                        ,range:: from_vector(v)
                );
    }

    struct {} ref_wraps;
    struct {} map;
    struct {} foreach;
    struct {} unzip_foreach;
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
    auto operator| (R r, decltype(map) ) {
        return temporary_tagged_holder<R, decltype(map)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(foreach) ) {
        return temporary_tagged_holder<R, decltype(foreach)> { std::move(r) };
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
    template<typename R, typename F>
    auto operator| (temporary_tagged_holder<R, decltype(unzip_map)> r_holder, F f)
    {
        return unzip_map_t<R,F>{ move(r_holder.m_r), std:: move(f) };
    }

    template<typename R, typename F>
    struct map_t {
        R m_r;
        F f;
        bool empty() const {
            return m_r.empty();
        }
        auto front_val() {
            return f( range:: front_val(m_r) );
        }
        void advance() { m_r.advance(); }
    };
    template<typename R, typename F>
    auto operator| (temporary_tagged_holder<R, decltype(map)> r_holder, F f)
    {
        return map_t<R,F>{ move(r_holder.m_r), std:: move(f) };
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
} // namespace range
#endif
