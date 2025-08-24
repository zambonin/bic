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
uintx inner_bin(const uint32_t n, const uint32_t k, const uint16_t d) {
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

uintx bin(const uint32_t n, const uint32_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(BIN_CACHE, GET_CACHE_BIN(n, k), inner_bin);
}

uintx alt_inner_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  intx rop = 0;
  intx inner = 0;
  uintx left = 0;
  uintx right = 0;

  uint32_t i = 1 + ((n + k - 1) / (d + 1));
  for (; i <= k; ++i) {
    left = inner_bin(k, i, d);
    right = inner_bin(i * (d + 1) - 1 - n, k - 1, d);
    inner = left * right;
    if ((k - i) & 1U) {
      inner = -inner;
    }
    rop += inner;
  }

  return (uintx)rop;
}

uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx inner = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = bin(k, i, d);
    right = bin(n - (d + 1) * i + k - 1, k - 1, d);
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
  return inner_bic_with_sums(n, k, d, NULL);
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
