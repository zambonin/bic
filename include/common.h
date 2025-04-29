#ifndef COMMON_H
#define COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

enum {
  BIT_LENGTH = 320,
  NS_TO_SEC = 1000000000,
};

#if defined(BOOST_FIX_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::checked_uint512_t;
using intx = boost::multiprecision::checked_int512_t;
#elif defined(BOOST_ARB_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::cpp_int;
using intx = boost::multiprecision::cpp_int;
#elif defined(BOOST_MPZ_INT)
#include <boost/multiprecision/gmp.hpp>

using uintx = boost::multiprecision::mpz_int;
using intx = boost::multiprecision::mpz_int;
#elif defined(BOOST_TOM_INT)
#include <boost/multiprecision/tommath.hpp>

using uintx = boost::multiprecision::tom_int;
using intx = boost::multiprecision::tom_int;
#endif

#if defined(BITINT)
typedef unsigned _BitInt(BIT_LENGTH) uintx;
typedef _BitInt(BIT_LENGTH) intx;
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/random.hpp>
#include <random>

static std::random_device rd;
#endif

typedef void (*unrank_func)(uint32_t *, const uint16_t, const uint16_t,
                            const uint16_t, const uintx);
typedef uintx (*rank_func)(const uint16_t, const uint16_t, const uint16_t,
                           const uint32_t *);
typedef void (*strategy_func)(const uint16_t, uint16_t *, const uint16_t,
                              uint16_t *);
typedef uintx (*math_func)(const uint16_t, const uint16_t, const uint16_t);

#endif
