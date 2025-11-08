#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stdlib.h>

typedef struct bic_ctx_s *bic_ctx_t;

#if defined(BOOST_FIX_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::checked_uint512_t;
using intx = boost::multiprecision::checked_int512_t;
static const double BIT_LENGTH = 512;
#elif defined(BOOST_ARB_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::cpp_int;
using intx = boost::multiprecision::cpp_int;
static const double BIT_LENGTH = INFINITY;
#elif defined(BOOST_MPZ_INT)
#include <boost/multiprecision/gmp.hpp>

using uintx = boost::multiprecision::mpz_int;
using intx = boost::multiprecision::mpz_int;
static const double BIT_LENGTH = INFINITY;
#elif defined(BOOST_TOM_INT)
#include <boost/multiprecision/tommath.hpp>

using uintx = boost::multiprecision::tom_int;
using intx = boost::multiprecision::tom_int;
static const double BIT_LENGTH = INFINITY;
#endif

#if defined(BITINT)
typedef unsigned _BitInt(BITINT) uintx;
typedef _BitInt(BITINT) intx;
static const double BIT_LENGTH = BITINT;
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

#if defined(BOOST_FIX_INT) || defined(BOOST_ARB_INT) ||                        \
    defined(BOOST_MPZ_INT) || defined(BOOST_TOM_INT)
#define uintx_alloc(count) (new uintx[count])
#define uintx_free(ptr) (delete[] ptr)
#define intx_alloc(count) (new intx[count])
#define intx_free(ptr) (delete[] ptr)
#elif defined(BITINT)
#define uintx_alloc(count) ((uintx *)calloc(count, sizeof(uintx)))
#define uintx_free(ptr) (free(ptr))
#define intx_alloc(count) ((intx *)calloc(count, sizeof(intx)))
#define intx_free(ptr) (free(ptr))
#endif

#endif
