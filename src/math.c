#include "math.h"
#include "cache.h"
#include "utils.h"

#if defined(BITINT)
long double lg(const uintx u) {
  // actually a truncated logarithm
  if (u == 0) {
    return -1.0;
  }
  for (uint16_t i = BIT_LENGTH - 1; i > 0; i--) {
    if ((u >> i) & 1) {
      return i;
    }
  }
  return 0;
}
#else
long double lg(const uintx u) {
  return (long double)log2(boost::multiprecision::cpp_bin_float_100(u));
}
#endif

uint16_t bits_fit_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return 1 + ((uint16_t)lg(inner_bic_with_sums(n, k, d, NULL, inner_bin)));
}

double asqrt(double x) {
  if (x < 0) {
    return -1.0;
  }
  if (x == 0) {
    return 0;
  }
  const double e = 0.000001;

  double z = x / 2.0;
  double y = 0.0;
  while (z - y > e) {
    y = z;
    z = (z + x / z) / 2.0;
  }
  return z;
}

// from FXT: aux0/binomial.h
uintx inner_bin(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  if (k > n) {
    return 0;
  }

  if ((k == 0) || (k == n)) {
    return 1;
  }

  uint16_t kk = k;
  if (2 * k > n) {
    kk = n - k;
  }

  uintx b = n - kk + 1;
  uintx f = b;

  for (uintx j = 2; j <= kk; ++j) {
    ++f;
    b *= f;
    b /= j;
  }

  return b;
}

uintx bin(const uint16_t n, const uint16_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(BIN_CACHE, GET_CACHE_BIN(n, k), inner_bin);
}

uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums, math_func bin_impl) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx inner = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = bin_impl(k, i, d);
    right = bin_impl(n - (d + 1) * i + k - 1, k - 1, d);
    inner = left * right;
    if (i & 1U) {
      inner = -inner;
    }
    rop += inner;

    if (partial_sums != NULL) {
      partial_sums[i] = inner;
    }
  }

  return (uintx)rop;
}

uintx inner_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return inner_bic_with_sums(n, k, d, NULL, bin);
}

uintx bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  if (cache_type == SMALL_COMB_CACHE) {
    uintx *row = GET_CACHE_SCOMB(0, k - 1);
    uint16_t left = (uint16_t)row[0];
    uint16_t right = (uint16_t)row[1];

    if (n < left || n > right) {
      return inner_bic(n, k, d);
    }
    return row[n - left + 2];
  }

  GET_CACHE_OR_CALC(COMB_CACHE, GET_CACHE_COMB(n, k), inner_bic);
}

uintx *inner_acc(const uint16_t n, const uint16_t k, const uint16_t d) {
  size_t length = d + 3;
  size_t i = 0;
  uintx *rop = (uintx *)calloc(length, sizeof(uintx));
  acc_cache_t.total_size += length * sizeof(uintx);
  uintx sum = 0;

  rop[0] = 0;
  for (; i <= min(n, d); ++i) {
    sum += bic(n - i, k, d);
    rop[i + 1] = sum;
  }
  rop[length - 1] = i;

  return rop;
}

uintx *acc(const uint16_t n, const uint16_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(ACC_COMB_CACHE, GET_CACHE_ACC(n, k), inner_acc);
}

uintx bic_acc(const uint16_t n, const uint16_t k, const uint16_t d,
              const uint16_t l) {
  if (cache_type == ACC_COMB_CACHE) {
    return acc(n, k, d)[l];
  }

  uint16_t j = min(k, n / (d + 1));
  uint16_t u;

  intx tmp = 0;
  intx rop = 0;

  for (uint16_t i = 0; i <= j; ++i) {
    u = n - (d + 1) * i + k;
    tmp = bin(k, i, d) * (bin(u, k, d) - bin(max(0, u - l), k, d));
    if (i & 1U) {
      tmp = -tmp;
    }
    rop += tmp;
  }

  return (uintx)rop;
}

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d) {
  uint16_t len = bits_fit_bic(n, k, d) / sizeof(uint64_t);
  len += (len == 0);

  uint8_t *message = (uint8_t *)calloc(len, sizeof(uint8_t));
  for (uint16_t i = 0; i < len; ++i) {
    message[i] = random();
  }
  uintx rank = 0;

#if defined(BITINT)
  unsigned char *ptr = (unsigned char *)&rank;
  for (uint16_t i = 0; i < len; ++i) {
    ptr[len - 1 - i] = message[i];
  }
#elif defined(BOOST_FIX_INT) || defined(BOOST_ARB_INT)
  boost::multiprecision::detail::import_bits_fast(rank, message, message + len);
#elif defined(BOOST_MPZ_INT)
  mpz_import(rank.backend().data(), len, 1, sizeof(uint8_t), 0, 0, message);
#elif defined(BOOST_TOM_INT)
  mp_err err =
      mp_unpack(&rank.backend().data(), len, 1, sizeof(uint8_t), 0, 0, message);
  (void)err;
#endif

  free(message);
  return rank % inner_bic(n, k, d);
}
