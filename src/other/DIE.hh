#ifdef __EMSCRIPTEN__
#define error_print_and_exit(...) ([&]()->bool{ std::cout << '\n' << "Exiting: "  __FILE__  ":" << __LINE__ << "\tError: " << __VA_ARGS__; std::cout << '\n'; exit(1);}())
#else
#define error_print_and_exit(...) ([&]()->bool{ std::cerr << '\n' << "Exiting: "  __FILE__  ":" << __LINE__ << "\tError: " << __VA_ARGS__; std::cerr << '\n'; exit(1);}())
#endif
#define DIE(...) error_print_and_exit(__VA_ARGS__)

#define WARNING(...) ([&]()->bool{ std::cerr << '\n' << "WARNING: " << __VA_ARGS__; std::cerr << '\n'; return false; }())
