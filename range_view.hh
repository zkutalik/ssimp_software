#ifndef AMD_RANGE_VIEW_HH
#define AMD_RANGE_VIEW_HH

#include"range.hh"
#include<functional>

namespace range {
namespace view {
    template<typename V>
    auto enumerate_vector(V && v) {
        return zip_val  (range:: ints       ( AMD_FORWARD(v).size())
                        ,range:: from_vector( AMD_FORWARD(v) )
                );
    }
    template<typename R>
    auto enumerate(R && r) {
        return zip_val  (range:: ints       (                )
                        ,AMD_FORWARD(r)
                );
    }

    struct {} ref_wraps;
    struct {} which;        // given a range of bools, return the indices (int64_t) of those that are true
    struct {} map;
    struct {} take;
    struct {} unzip_filter;
    struct {} unzip_map;

    template<typename R>
    auto operator| (R r, decltype(ref_wraps) );
    template<typename R>
    auto operator| (R r, decltype(which) );

    template<typename R>
    struct ref_wraps_impl {
        R m_range_with_front_ref;

        bool            empty()         const   { return   range:: empty    (m_range_with_front_ref)  ; }
        auto            front_val()     const   { return std:: ref( range:: front_ref(m_range_with_front_ref) ); }
        void            advance()               { range:: advance(m_range_with_front_ref); }
    };

    template<typename R>
    auto operator| (R r, decltype(ref_wraps) ) {
        return ref_wraps_impl<R> {std:: move(r)};
    }

    namespace impl {
        static
        constexpr int64_t special_index_notevenstarted = -1;
        static
        constexpr int64_t special_index_empty = -1337;
        template<typename R>
        struct which_impl {
            mutable
            R m_r;
            mutable
            int64_t     m_index; /* -1 means we haven't done anything yet
                                  * -1337 means we're at empty
                                  *  anything else means  true==front_val(m_r)
                                  */

            bool            empty()         const   {
                if(m_index == special_index_notevenstarted) {
                    real_advance();
                    return empty();
                }
                assert(m_index != special_index_notevenstarted);
                return m_index == special_index_empty;
            }
            uint64_t        front_val()     const   {
                if(m_index == special_index_notevenstarted) {
                    real_advance();
                    return this->front_val();
                }
                assert(m_index != special_index_notevenstarted);
                assert(m_index != special_index_empty);
                return m_index;
            }
            void            real_advance()     const {
                assert(m_index != special_index_empty);


                bool b;
                do{
                    if( range::empty(m_r) ) {
                        m_index = special_index_empty;
                        return;
                    }

                    b = range:: pull(m_r);
                    ++m_index;
                    // that pre-increment above only works because special_index_notevenstarted==-1
                    static_assert(special_index_notevenstarted == -1 ,"");
                } while(!b);
            }
            void    advance() {
                real_advance();
            }
        };
    }
    template<typename R>
    auto operator| (R r, decltype(which) ) {
        return impl:: which_impl<R> {std:: move(r), impl::special_index_notevenstarted};
    }

    template<typename R>
    auto operator| (R r, decltype(unzip_filter) ) {
        return detail:: temporary_tagged_holder<R, decltype(unzip_filter)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(unzip_map) ) {
        return detail:: temporary_tagged_holder<R, decltype(unzip_map)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(map) ) {
        return detail:: temporary_tagged_holder<R, decltype(map)> { std::move(r) };
    }
    template<typename R>
    auto operator| (R r, decltype(take) ) {
        return detail:: temporary_tagged_holder<R, decltype(take)> { std::move(r) };
    }

    template<typename R, typename F>
    struct unzip_filter_t {
        R m_r;
        F m_f;

        static_assert(!std::is_reference<R>{} ,"");
        static_assert(!std::is_reference<F>{} ,"");

        unzip_filter_t(R r, F f) : m_r(std::move(r)), m_f(std::move(f))
        {
            skip_filter_mismatches();
        }
        bool empty() const {
            return m_r.empty();
        }
        auto front_val() -> AMD_RANGE_DECLTYPE_AND_RETURN( range:: front_val(m_r) )
        void advance() {
            m_r.advance();
            skip_filter_mismatches();
        }
        void skip_filter_mismatches() {
            while   (   !empty()
                     && !utils:: apply  ( std::forward<F>(m_f)
                                        , range:: front_val(m_r)
                                        )) {
                m_r.advance();
            }
        }
    };
    template<typename R, typename F>
    auto operator| (detail:: temporary_tagged_holder<R, decltype(unzip_filter)> r_holder, F f)
    {
        return unzip_filter_t<R,F>{ move(r_holder.m_r), std:: move(f) };
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
    auto operator| (detail:: temporary_tagged_holder<R, decltype(unzip_map)> r_holder, F f)
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
        auto front_val() const {
            return f( range:: front_val(m_r) );
        }
        void advance() { m_r.advance(); }
    };
    template<typename R, typename F>
    auto operator| (detail:: temporary_tagged_holder<R, decltype(map)> r_holder, F f)
    {
        return map_t<R,F>{ std:: move(r_holder.m_r), std:: move(f) };
    }

    template<typename R>
    struct take_t {
        R m_r;
        int64_t     m_how_many;
        bool empty() const {
            return m_how_many==0 || m_r.empty();
        }
        auto front_val() const {
            return range:: front_val(m_r);
        }
        void advance() { --m_how_many; m_r.advance(); }
    };
    template<typename R>
    auto operator| (detail:: temporary_tagged_holder<R, decltype(take)> r_holder, int64_t how_many)
    {
        return take_t<R>{ std:: move(r_holder.m_r), how_many };
    }

} // namespace view
} // namespace range
#endif
