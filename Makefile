# required libraries:
# Boost >= 1.56
# Boost.SIMD

# define library paths here:
# BOOST_ROOT = /path/to/boost
# SIMD_INCLUDE = /path/to/boost-simd/include
# odeint is part of the latest boost releases
# ODEINT_ROOT = ${BOOST_ROOT}

ARCH = Host
# ARCH = AVX

CXX = icpc
CC = icpc
CXXFLAGS = -x${ARCH} -Ofast -fno-alias -ip -inline-forceinline -DNDEBUG -std=c++0x -I${ODEINT_ROOT} -I${BOOST_ROOT} -I${SIMD_INCLUDE}

#g++ compiler

# CXX = g++	
# CC = g++
# TUNE = -mtune=native -march=native
# # TUNE = -mtune=btver1 -march=btver1
# # TUNE = -march=corei7-avx
# # TUNE = -mtune=sse42 

# CXXFLAGS = -O3 -ffast-math -std=c++0x ${TUNE} -DNDEBUG -I${ODEINT_ROOT} -I${BOOST_ROOT} -I${SIMD_INCLUDE}

all: roessler_std roessler_cache_opt roessler_cache_opt_simd
