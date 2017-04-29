#ifndef AMD_RANGE_ACTION_HH
#define AMD_RANGE_ACTION_HH

#include"range.hh"

#include <type_traits>
#include <vector>

namespace range {
namespace action {

    struct {} collect;
    struct {} unzip_collect_transpose;

    template<typename R>
    auto operator| (R r, decltype(collect) )
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
    {
        return std::make_tuple( std:: get<Is>( move( r.m_ranges)) | collect ... );
    }
    template<typename R>
    auto operator| (R r, decltype(unzip_collect_transpose) )
    {
        return operator_pipe_impl( move(r), unzip_collect_transpose
                , std:: make_index_sequence< r.width_v >{}
                );
    }
}
}

#endif
