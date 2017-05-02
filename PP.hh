#ifndef MACROS_HPP__
#define MACROS_HPP__

#define PP1(x)                 std :: cout << #x << ":" << (x) << std :: endl
#define PP2(x,y)               std :: cout << #x << ',' << #y                                           << ":\t" << (x) << " , " << (y) << std :: endl
#define PP3(x,y,z)             std :: cout << #x << ',' << #y << ',' << #z                              << ":\t" << (x) << " , " << (y) << " , " << (z) << std :: endl
#define PP4(x,y,z,w)           std :: cout << #x << ',' << #y << ',' << #z << ',' << #w                 << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << std :: endl
#define PP5(x,y,z,w,v)         std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v    << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << std :: endl
#define PP6(x,y,z,w,v,u)       std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u) << std :: endl
#define PP7(x,y,z,w,v,u,t)     std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u << ',' << #t      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u) << " , " << (t) << std :: endl
#define PP8(x,y,z,w,v,u,t,s)   std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u << ',' << #t << ',' << #s      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u) << " , " << (t) << " , " << (s) << std :: endl
#define PP9(x,y,z,w,v,u,t,s,r)   std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u <<','  <<#t  <<','  <<#s  <<','  <<#r      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u)             <<" , "<<(t) <<" , "<<(s) <<" , "<<(r) << std :: endl
#define PP10(x,y,z,w,v,u,t,s,r,q) do {\
	std::cout<< #x  << ',' <<#y  << ',' <<#z  << ',' <<#w  << ',' <<#v  << ',' <<#u  << ',' <<#t  << ',' <<#s  << ',' <<#r  << ',' <<#q      \
	<< ":\t" << (x) <<" , "<<(y) <<" , "<<(z) <<" , "<<(w) <<" , "<<(v) <<" , "<<(u) <<" , "<<(t) <<" , "<<(s) <<" , "<<(r) <<" , "<<(q) <<std::endl; \
} while(0)


#define GET_ARG_11(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,N,...) N
#define COUNT_MACRO_ARGS(...) GET_ARG_11( __VA_ARGS__, 10,9,8,7,6,5,4,3,2,1,I_CANNOT_SEE_ZERO_ARGS)

#define SELECT_PP_IMPL( n ) PP ## n
#define SELECT_PP( n ) SELECT_PP_IMPL(n)
#define PP(...) do{ SELECT_PP( COUNT_MACRO_ARGS(__VA_ARGS__) )(__VA_ARGS__) ;}while(0)

#define PPL(...) do{ cout << __FILE__ ":" << __LINE__ << "\t   "; PP(__VA_ARGS__); } while(0)

#endif
