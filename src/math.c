#include "math.h"
#include "cache.h"
#include "primes.h"

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

uint64_t kummer(const uint64_t n, const uint64_t k, const uint64_t p) {
  uint64_t n_ = n;
  uint64_t k_ = k;

  uint64_t c = 0;
  uint64_t r = 0;
  while (n_ >= p) {
    if ((n_ % p) < (k_ % p) + r) {
      r = 1;
      c++;
    } else {
      r = 0;
    }

    n_ /= p;
    k_ /= p;
  }

  return c + (n_ < k_ + r);
}

uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  (void)ctx;

  if (k > n) {
    return 0;
  }
  if (k == 0 || k == n) {
    return 1;
  }

  const uint64_t kb = (k > n / 2) ? n - k : k;
  const uint64_t nk = n - kb;

  uintx res = 1;
  for (size_t i = 0; i < N_PRIMES; i++) {
    const uint64_t p = PRIMES[i];
    if (p > n) {
      break;
    }

    if (p > nk) {
      res *= (uintx)p;
    } else if (p > n / 2) {
      continue;
    } else if (p * p > n) {
      if ((n % p) < (kb % p)) {
        res *= (uintx)p;
      }
    } else {
      const uint64_t e = kummer(n, kb, p);
      if (e) {
        res *= ipow(p, e);
      }
    }
  }

  return res;
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

  return (rop < 0) ? 0 : (uintx)rop;
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

  return (rop < 0) ? 0 : (uintx)rop;
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

  return (rop < 0) ? 0 : (uintx)rop;
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

  uintx total_sum = 0;
  for (uint32_t i = 0; i <= l; ++i) {
    if (n < i) {
      break;
    }

    uintx term = ctx->comp(n - i, k, d, ctx);
    if (p > 0) {
      term *= ipow(n - i, p);
    }
    total_sum += term;
  }

  return total_sum;
}

void S0_to_S2(uintx *S0, uintx *S1, uintx *S2, const uint32_t n,
              const uint32_t k, const uint32_t d, const uint32_t l,
              const bic_ctx_t ctx) {
  if (k == 0) {
    *S0 = *S1 = *S2 = (n == 0 && l > 0) ? 1 : 0;
    return;
  }

  *S0 = *S1 = *S2 = 0;
  for (uint32_t i = 0; i <= l; ++i) {
    if (n < i) {
      break;
    }
    uintx term = ctx->comp(n - i, k, d, ctx);
    if (term > 0) {
      *S0 += term;
      *S1 += term * (n - i);
      *S2 += term * ipow(n - i, 2);
    }
  }
}

double exp_total_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                     const uint16_t l, const bic_ctx_t ctx) {
  if (l == 0) {
    return (double)n;
  }

  uintx S0 = 0;
  uintx S1 = 0;
  uintx S2 = 0;

  S0_to_S2(&S0, &S1, &S2, n, k, d, l, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  return (double)S1 / (double)S0;
}

double var_total_sum(const uint16_t n, const uint16_t k, const uint16_t d,
                     const uint16_t l, const bic_ctx_t ctx) {
  if (l == 0) {
    return 0.0;
  }

  uintx S0 = 0;
  uintx S1 = 0;
  uintx S2 = 0;

  S0_to_S2(&S0, &S1, &S2, n, k, d, l, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  double exp_val = (double)S1 / (double)S0;
  return dmax(0.0, ((double)S2 / (double)S0) - (exp_val * exp_val));
}

double var_single_part(const uint16_t n, const uint16_t k, const uint16_t d,
                       const uint16_t l, const bic_ctx_t ctx) {
  if (k == 0) {
    return 0.0;
  }

  uintx S0 = S_p(n, k, d, l, 0, ctx);
  if (S0 == 0) {
    return 0.0;
  }

  double lhs = 0.0;
  intx prev_S0 = (intx)S_p(n, k - 1, d, l, 0, ctx);
  for (uint32_t y = 0; y <= min(n, d); ++y) {
    lhs += ((double)y * y) * (double)prev_S0;
    if (y == min(n, d)) {
      break;
    }
    const intx ny = (intx)ctx->comp(n - y, k - 1, d, ctx);
    const int32_t term = n - y - l - 1;
    intx nyl = 0;
    if (term >= 0) {
      nyl = (intx)ctx->comp(term, k - 1, d, ctx);
    }
    prev_S0 = prev_S0 - ny + nyl;
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
