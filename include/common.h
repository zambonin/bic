#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

static const uint32_t NS_TO_SEC = 1000000000;

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

typedef struct {
  void (*unrank)(uint32_t *, const uint16_t, const uint16_t, const uint16_t,
                 const uintx);
  uintx (*rank)(const uint16_t, const uint16_t, const uint16_t,
                const uint32_t *);
} order;

typedef void (*strategy_func)(const uint16_t, uint16_t *, const uint16_t,
                              uint16_t *);

#endif
