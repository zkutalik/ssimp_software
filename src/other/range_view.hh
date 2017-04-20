#include"range.hh"
#include<functional>

namespace view {
    struct {} ref_wraps;
    struct {} foreach;
    struct {} unzip_foreach;
    struct {} collect;

    template<typename R, typename tag>
    struct temporary_tagged_holder {
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
    void operator| (temporary_tagged_holder<R, decltype(foreach)> r_holder, F && f)
    {
        while(!r_holder.m_r.empty()) {
            std::forward<F>(f)(
                range:: pull(r_holder.m_r)
            );
        }
    }

} // namespace view
