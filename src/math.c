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
  while (z - y > e) {
    y = z;
    z = (z + x / z) / 2.0;
  }
  return z;
}

uintx compute_bin(const bic_ctx_t *ctx, const uint32_t n, const uint32_t k) {
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

uintx alt_compute_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint32_t i = 1 + ((n + k - 1) / (d + 1));
  for (; i <= k; ++i) {
    left = compute_bin(NULL, k, i);
    right = compute_bin(NULL, i * (d + 1) - 1 - n, k - 1);
    compute = left * right;
    if ((k - i) & 1U) {
      compute = -compute;
    }
    rop += compute;
  }

  return (uintx)rop;
}

uintx compute_bic_with_sums(const bic_ctx_t *ctx, const uint16_t n,
                            const uint16_t k, const uint16_t d,
                            intx *partial_sums) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = ctx->bin(ctx, k, i);
    right = ctx->bin(ctx, n - (d + 1) * i + k - 1, k - 1);
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

uintx compute_bic(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                  const uint16_t d) {
  return compute_bic_with_sums(ctx, n, k, d, NULL);
}

void compute_acc(uintx *rop, const bic_ctx_t *ctx, const uint16_t n,
                 const uint16_t k, const uint16_t d) {
  const size_t length = d + 3;
  uint16_t i = 0;
  uintx sum = 0;

  rop[0] = 0;
  for (; i <= min(n, d); ++i) {
    sum += ctx->comp(ctx, n - i, k, d);
    rop[i + 1] = sum;
  }
  rop[length - 1] = i;
}

uintx compute_dir(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                  const uint16_t d, const uint16_t l) {
  uint16_t j = min(k, n / (d + 1));
  uint16_t u;

  intx tmp = 0;
  intx rop = 0;

  for (uint16_t i = 0; i <= j; ++i) {
    u = n - (d + 1) * i + k;
    tmp = ctx->bin(ctx, k, i) *
          (ctx->bin(ctx, u, k) - ctx->bin(ctx, max(0, u - l), k));
    if (i & 1U) {
      tmp = -tmp;
    }
    rop += tmp;
  }

  return (uintx)rop;
}

uintx scomb(const uint32_t n, const uint32_t k) {
  return compute_bin(NULL, n, k);
}

uintx power(uintx base, int exp) {
  uintx res = 1;
  for (int i = 0; i < exp; ++i) {
    res *= base;
  }
  return res;
}

uintx W(int s, int c, int y, int l) {
  if (s < 0 || c < 0 || y < 0) {
    return 0;
  }
  if (s == 0) {
    return scomb(y + 1, c + 1) - scomb(y - l, c + 1);
  }
  return ((uintx)c * W(s - 1, c, y, l)) +
         ((uintx)(c + 1) * W(s - 1, c + 1, y, l));
}

uintx S_p(uint32_t n, uint32_t k, uint32_t d, uint32_t l, uint32_t p) {
  if (k == 0) {
    return (n == 0 && l > 0) ? 1 : 0;
  }

  uintx total_sum = 0;
  uint32_t j = (d > 0) ? min(k, n / (d + 1)) : k;

  for (uint32_t i = 0; i <= j; ++i) {
    uintx u_i = (uintx)(d + 1) * i - k + 1;
    uintx inner = 0;
    for (uint32_t s = 0; s <= p; ++s) {
      inner += scomb(p, s) * power(u_i, p - s) * W(s, k - 1, n - u_i, l);
    }

    uintx term = scomb(k, i) * inner;
    if (i % 2 != 0) {
      total_sum -= term;
    } else {
      total_sum += term;
    }
  }

  return total_sum;
}

uintx _S0(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  return S_p(n, k, d, l, 0);
}

uintx _S1(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  return S_p(n, k, d, l, 1);
}

uintx _S2(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  return S_p(n, k, d, l, 2);
}

double exp_total_sum(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  if (l == 0) {
    return (double)n;
  }
  uintx S0 = _S0(n, k, d, l);
  if (S0 == 0) {
    return 0.0;
  }
  return (double)_S1(n, k, d, l) / (double)S0;
}

double var_total_sum(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  if (l == 0) {
    return 0.0;
  }
  uintx S0 = _S0(n, k, d, l);
  if (S0 == 0) {
    return 0.0;
  }
  double exp_val = exp_total_sum(n, k, d, l);
  return dmax(0.0,
              ((double)_S2(n, k, d, l) / (double)S0) - (exp_val * exp_val));
}

double var_single_part(uint32_t n, uint32_t k, uint32_t d, uint32_t l) {
  if (k == 0) {
    return 0.0;
  }
  uintx S0 = _S0(n, k, d, l);
  if (S0 == 0) {
    return 0.0;
  }

  double lhs = 0.0;
  for (uint32_t y = 0; y <= min(n, d); ++y) {
    lhs += ((double)y * y) * (double)S_p(n - y, k - 1, d, l, 0);
  }
  lhs /= (double)S0;

  double exp_val = exp_total_sum(n, k, d, l) / k;
  double rhs = exp_val * exp_val;

  return dmax(0.0, lhs - rhs);
}

double exp_part_sum(uint32_t n, uint32_t k, uint32_t d, uint32_t j,
                    uint32_t l) {
  if (k <= 0) {
    return 0.0;
  }
  return (double)j * exp_total_sum(n, k, d, l) / k;
}

double var_part_sum(uint32_t n, uint32_t k, uint32_t d, uint32_t j,
                    uint32_t l) {
  if (k <= 1) {
    return var_total_sum(n, k, d, l);
  }

  double lhs = ((double)j * (k - j) / (k - 1)) * var_single_part(n, k, d, l);
  double rhs =
      ((double)j * (j - 1) / ((double)k * (k - 1))) * var_total_sum(n, k, d, l);

  return dmax(0.0, lhs + rhs);
}

double stddev_part_sum(uint32_t n, uint32_t k, uint32_t d, uint32_t j,
                       uint32_t l) {
  return asqrt(var_part_sum(n, k, d, j, l));
}
