#ifdef __GNUC__
	#define assert(condition) do { if(condition) {} else { stop(std::string("assert() failure: (" #condition ") in ") + __PRETTY_FUNCTION__); } }while(0)
#else
	#define assert(condition) do { if(condition) {} else { stop(            "assert() failure: (" #condition ")"                           ); } }while(0)
#endif
