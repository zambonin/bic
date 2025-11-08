#include "math.h"
#include "cache.h"

uint32_t min(const uint32_t a, const uint32_t b) { return (a < b) ? a : b; }

int32_t max(const int32_t a, const int32_t b) { return (a > b) ? a : b; }

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

uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  (void)ctx;
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

uintx compute_bic_no_cache(const uint16_t n, const uint16_t k,
                           const uint16_t d) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = compute_bin(k, i, NULL);
    right = compute_bin(n - (d + 1) * i + k - 1, k - 1, NULL);
    compute = left * right;
    if (i & 1U) {
      compute = -compute;
    }
    rop += compute;
  }

  return (uintx)rop;
}

uintx compute_bic_with_sums(const uint16_t n, const uint16_t k,
                            const uint16_t d, intx *partial_sums,
                            const bic_ctx_t ctx) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = ctx->bin(k, i, ctx);
    right = ctx->bin(n - (d + 1) * i + k - 1, k - 1, ctx);
    compute = left * right;
    if (i & 1U) {
      compute = -compute;
    }
    rop += compute;

    if (partial_sums != NULL) {
      partial_sums[i] = compute;
    }
  }

  return (uintx)rop;
}

uintx compute_bic(const uint16_t n, const uint16_t k, const uint16_t d,
                  const bic_ctx_t ctx) {
  return compute_bic_with_sums(n, k, d, NULL, ctx);
}

uint16_t compute_acc(uintx *rop, const uint16_t n, const uint16_t k,
                     const uint16_t d, const bic_ctx_t ctx) {
  uint16_t i = 0;
  uintx sum = 0;

  rop[0] = 0;
  for (; i <= min(n, d); ++i) {
    sum += ctx->comp(n - i, k, d, ctx);
    rop[i + 1] = sum;
  }

  return i;
}

uintx compute_dir(const uint16_t n, const uint16_t k, const uint16_t d,
                  const uint16_t l, const bic_ctx_t ctx) {
  uint16_t j = min(k, n / (d + 1));
  uint16_t u;

  intx tmp = 0;
  intx rop = 0;

  for (uint16_t i = 0; i <= j; ++i) {
    u = n - (d + 1) * i + k;
    tmp = ctx->bin(k, i, ctx) *
          (ctx->bin(u, k, ctx) - ctx->bin(max(0, u - l), k, ctx));
    if (i & 1U) {
      tmp = -tmp;
    }
    rop += tmp;
  }

  return (uintx)rop;
}
