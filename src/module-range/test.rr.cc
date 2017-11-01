#include "rr.hh"
#include "../bits.and.pieces/PP.hh"
#include "../bits.and.pieces/utils.hh"
#include<iostream>
#include<vector>
using std:: vector;
using std:: string;
using namespace rr;
using utils:: operator<<;
int main () {
    auto r_i = ints(3);
    while(!empty(r_i)) {
        PP(front_val(r_i));
        advance(r_i);
    }

    r_i = ints(4);
    for(auto i : r_i) { PP(i); }

    vector<string> v{"hi", "world", "of", "ranges"};
    auto v_r = as_range(v);
    static_assert(                                      rr:: is_range_v< decltype(v_r) >    , "");
    static_assert(                                    ! rr:: is_range_v< decltype(v  ) >    , "");
    while(!empty(v_r)) {
        PP(front_val(v_r));
        advance(v_r);
    }
    auto mapped = v |map_range|[](auto && x)->int{
        return x.length();
    };
    while(!empty(mapped)) {
        PP(front_val(mapped));
        advance(mapped);
    }

    auto o = v |map_range|[](auto && x)->int{
        return -x.length();
    };
    for(;!empty(o);advance(o)) {
        PP(front_val(o));
    }
    auto collected = v |map_collect|[](auto && x) { return 0.5+x.length(); };
    auto recollected = collected |map_range|[](auto && x){return -x;} |collect;
    PP(collected, recollected);

}
