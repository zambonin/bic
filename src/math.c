#include "math.h"
#include "cache.h"

uint32_t min(const uint32_t a, const uint32_t b) { return (a < b) ? a : b; }

int32_t max(const int32_t a, const int32_t b) { return (a > b) ? a : b; }

long double dmax(const long double a, const long double b) {
  return (a > b) ? a : b;
}

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
  while (y == 0 || (y - z) > e) {
    y = z;
    z = (z + x / z) / 2.0;
  }
  return z;
}

uintx ipow(const uint64_t b, const uint64_t e) {
  uint64_t base = b;
  uint64_t exp = e;

  uintx result = 1;
  while (exp) {
    if (exp & 1) {
      result *= base;
    }

    base *= base;
    exp >>= 1;
  }

  return result;
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

uintx W(const int16_t s, const int16_t c, const int16_t y, const int16_t l,
        const bic_ctx_t ctx) {
  if (s < 0 || c < 0 || y < 0) {
    return 0;
  }

  if (s == 0) {
    return (uintx)(ctx->bin(y + 1, c + 1, ctx) - ctx->bin(y - l, c + 1, ctx));
  }

  uintx sum = 0;
  for (int i = max(c, y - l); i <= y; ++i) {
    uintx pow_i_s = ipow(i, s);
    uintx bin_i_c = ctx->bin(i, c, ctx);
    sum += pow_i_s * bin_i_c;
  }
  return sum;
}

uintx S_p(const uint16_t n, const uint16_t k, const uint16_t d,
          const uint16_t l, const uint16_t p, const bic_ctx_t ctx) {
  if (k == 0) {
    return (n == 0 && l > 0) ? 1 : 0;
  }

  intx total_sum = 0;
  uint32_t j = (d > 0) ? min(k, n / (d + 1)) : k;

  for (uint32_t i = 0; i <= j; ++i) {
    int32_t u_i = (d + 1) * i - k + 1;
    intx inner = 0;
    for (uint32_t s = 0; s <= p; ++s) {
      inner += (intx)ctx->bin(p, s, ctx) * ipow(u_i, p - s) *
               W(s, k - 1, n - u_i, l, ctx);
    }

    intx term = (intx)ctx->bin(k, i, ctx) * inner;
    if (i % 2 != 0) {
      total_sum -= term;
    } else {
      total_sum += term;
    }
  }

  return (uintx)total_sum;
}

uintx _S0(const uint16_t n, const uint16_t k, const uint16_t d,
          const uint16_t l, const bic_ctx_t ctx) {
  return S_p(n, k, d, l, 0, ctx);
}

uintx _S1(const uint16_t n, const uint16_t k, const uint16_t d,
          const uint16_t l, const bic_ctx_t ctx) {
  return S_p(n, k, d, l, 1, ctx);
}

uintx _S2(const uint16_t n, const uint16_t k, const uint16_t d,
          const uint16_t l, const bic_ctx_t ctx) {
  return S_p(n, k, d, l, 2, ctx);
}

double exp_total_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                     const uint16_t l, bic_ctx_t ctx) {
  if (l == 0) {
    return (double)n;
  }

  uintx S0 = _S0(n, k, d, l, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  return (double)_S1(n, k, d, l, ctx) / (double)S0;
}

double var_total_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                     const uint16_t l, bic_ctx_t ctx) {
  if (l == 0) {
    return 0.0;
  }

  uintx S0 = _S0(n, k, d, l, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  double exp_val = exp_total_sum(n, k, d, l, ctx);
  return dmax(0.0, ((double)_S2(n, k, d, l, ctx) / (double)S0) -
                       (exp_val * exp_val));
}

double var_single_part(const uint16_t n, const uint16_t k, const uint16_t d,
                       const uint16_t l, const bic_ctx_t ctx) {
  if (k == 0) {
    return 0.0;
  }

  uintx S0 = _S0(n, k, d, l, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  double lhs = 0.0;
  for (uint32_t y = 0; y <= min(n, d); ++y) {
    lhs += ((double)y * y) * (double)S_p(n - y, k - 1, d, l, 0, ctx);
  }
  lhs /= (double)S0;

  double exp_val = exp_total_sum(n, k, d, l, ctx) / k;
  double rhs = exp_val * exp_val;

  return dmax(0.0, lhs - rhs);
}

double exp_part_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                    const uint16_t j, const uint16_t l, const bic_ctx_t ctx) {
  if (k <= 0) {
    return 0.0;
  }
  return (double)j * exp_total_sum(n, k, d, l, ctx) / k;
}

double var_part_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                    const uint16_t j, const uint16_t l, const bic_ctx_t ctx) {
  if (k <= 1) {
    return var_total_sum(n, k, d, l, ctx);
  }

  double lhs =
      ((double)j * (k - j) / (k - 1)) * var_single_part(n, k, d, l, ctx);
  double rhs = ((double)j * (j - 1) / ((double)k * (k - 1))) *
               var_total_sum(n, k, d, l, ctx);

  return dmax(0.0, lhs + rhs);
}

double stddev_part_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                       const uint16_t j, const uint16_t l,
                       const bic_ctx_t ctx) {
  return asqrt(var_part_sum(n, k, d, j, l, ctx));
}
