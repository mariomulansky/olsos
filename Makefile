
ARCH = Host
# ARCH = AVX

CXX = icpc
CC = icpc
CXXFLAGS = -x${ARCH} -Ofast -fno-alias -ip -inline-forceinline -DNDEBUG -std=c++0x -I${ODEINT_ROOT} -I${BOOST_ROOT} -I${SIMD_INCLUDE}

#g++ compiler

CXX = g++	
CC = g++
TUNE = -mtune=native -march=native
# TUNE = -mtune=btver1 -march=btver1
# TUNE = -march=corei7-avx
# TUNE = -mtune=sse42 

CXXFLAGS = -O3 -ffast-math -std=c++0x ${TUNE} -DNDEBUG -I${ODEINT_ROOT} -I${BOOST_ROOT} -I${SIMD_INCLUDE}
