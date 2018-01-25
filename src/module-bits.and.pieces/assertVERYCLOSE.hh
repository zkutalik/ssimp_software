#ifndef HH_assertVERYCLOSE_HH
#define HH_assertVERYCLOSE_HH

#define VERYCLOSE(a,b) (1e-06 > std::fabs((a)-(b)))
#define assertVERYCLOSE(a,b) assert(VERYCLOSE(a,b))

#endif
