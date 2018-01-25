#ifdef __GNUC__
	#define STRINGIFY_FOR_LINE_NUMBERS_(x) #x
	#define STRINGIFY_FOR_LINE_NUMBERS(x) STRINGIFY_FOR_LINE_NUMBERS_(x)
	#define assert(condition) do { if(condition) {} else { stop(std::string("assert() failure: (" #condition ") at [" __FILE__ ":" STRINGIFY_FOR_LINE_NUMBERS(__LINE__) "] in ") + __PRETTY_FUNCTION__); } }while(0)
#else
	#define assert(condition) do { if(condition) {} else { stop(            "assert() failure: (" #condition ")"                           ); } }while(0)
#endif
