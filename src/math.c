#include "math.h"
#include "primes.h"
#include <errno.h>
#include <float.h>

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
  for (uint32_t i = BIT_LENGTH - 1; i > 0; i--) {
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

double asqrt(const double x) {
  if (x < 0.0) {
    return -1.0;
  }
  if (x == 0.0 || x != x || x > DBL_MAX) {
    return x;
  }

  double y = x < 1.0 ? 1.0 : x;
  for (;;) {
    const double z = 0.5 * y + 0.5 * (x / y);
    if (z >= y) {
      return 0.5L * y + 0.5L * ((long double)x / y);
    }
    y = z;
  }
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

static uint32_t gcd32(uint32_t a, uint32_t b) {
  while (b != 0) {
    const uint32_t rem = a % b;
    a = b;
    b = rem;
  }
  return a;
}

#if defined(BITINT) || defined(BOOST_FIX_INT)
static uintx max_exact_value(void) {
  uintx limit = 1;
  limit <<= (uint32_t)BIT_LENGTH - 1;
  return limit - 1;
}

static bool binomial_bounded(const uint32_t n, const uint32_t k,
                             const uintx limit, uintx *out) {
  const uint32_t choose = min(k, n - k);
  uintx result = 1;
  for (uint32_t i = 1; i <= choose; ++i) {
    uint32_t numerator = n - choose + i;
    uint32_t denominator = i;
    const uint32_t common = gcd32(numerator, denominator);
    numerator /= common;
    denominator /= common;
    result /= denominator;
    if (result > limit / numerator) {
      return false;
    }
    result *= numerator;
  }
  if (out != NULL) {
    *out = result;
  }
  return true;
}
#endif

bool bic_params_supported(const uint32_t n, const uint32_t k, const uint32_t d,
                          const bic_cache_t cache) {
  if (cache >= BIC_CACHE_LENGTH || n > INT32_MAX || k > INT32_MAX ||
      d > INT32_MAX || (uint64_t)n + k + 1 > UINT32_MAX) {
    return false;
  }
  if (k == 0) {
    return n == 0;
  }

#if defined(BITINT) || defined(BOOST_FIX_INT)
  const uintx limit = max_exact_value();
  const uint32_t nk = n + k;
  const uint32_t terms = min(k, n / (d + 1));
  uintx choose_factor, bin_factor;
  if (!binomial_bounded(k, min(terms, k / 2), limit, &choose_factor) ||
      !binomial_bounded(nk, k, limit, &bin_factor) ||
      bin_factor > limit / choose_factor ||
      bin_factor * choose_factor > limit / (terms + 1)) {
    return false;
  }

  if (cache != BIC_CACHE_NONE &&
      !binomial_bounded(nk, min(k, nk / 2), limit, &bin_factor)) {
    return false;
  }
  if (cache >= BIC_CACHE_ACC) {
    if (!binomial_bounded(nk + 1, k + 1, limit, &bin_factor) ||
        bin_factor > limit / choose_factor ||
        bin_factor * choose_factor > limit / (terms + 1)) {
      return false;
    }
  }
#endif
  return true;
}

bool bic_dense_comb_params_supported(const uint32_t n, const uint32_t k,
                                     const uint32_t d) {
  if (n > INT32_MAX || k > INT32_MAX || d > INT32_MAX ||
      (uint64_t)n + k + 1 > UINT32_MAX) {
    return false;
  }
  if (k == 0) {
    return n == 0;
  }

#if defined(BITINT) || defined(BOOST_FIX_INT)
  const uintx limit = max_exact_value();
  const uint32_t base = d + 1;
  uintx bound = 1;
  for (uint32_t i = 1; i < k; ++i) {
    if (bound > limit / base) {
      return false;
    }
    bound *= base;
  }
#endif
  return true;
}

bool bic_cache_moments_supported(const uint32_t n, const uint32_t k,
                                 const uint32_t d) {
  if (!bic_dense_comb_params_supported(n, k, d)) {
    return false;
  }

#if defined(BITINT) || defined(BOOST_FIX_INT)
  const uintx limit = max_exact_value();
  const uint32_t sum = n == 0 ? 1 : n;
  const uint32_t base = d + 1;
  uintx bound = 1;
  for (uint32_t i = 0; i < 2; ++i) {
    if (bound > limit / sum) {
      return false;
    }
    bound *= sum;
  }
  for (uint32_t i = 0; i < k; ++i) {
    if (bound > limit / base) {
      return false;
    }
    bound *= base;
  }
#endif
  return true;
}

static uintx compute_bin_without_primes(const uint32_t n, const uint32_t k) {
  const uint32_t choose = min(k, n - k);
  uintx result = 1;

  for (uint32_t i = 1; i <= choose; ++i) {
    uint32_t numerator = n - choose + i;
    uint32_t denominator = i;
    const uint32_t common = gcd32(numerator, denominator);
    numerator /= common;
    denominator /= common;

    result /= denominator;
    result *= numerator;
  }

  return result;
}

// Goetgheluck (1987): doi:10.1080/00029890.1987.12000648
uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  (void)ctx;

  if (k > n) {
    return 0;
  }
  if (k == 0 || k == n) {
    return 1;
  }

#if defined(BITINT) || defined(BOOST_FIX_INT)
  if (n >= BIT_LENGTH && !binomial_bounded(n, k, max_exact_value(), NULL)) {
    errno = ERANGE;
    return 0;
  }
#endif

  if (n > PRIMES[N_PRIMES - 1]) {
    return compute_bin_without_primes(n, k);
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
        const uintx power = ipow(p, e);
        res *= power;
      }
    }
  }

  return res;
}

static uint32_t reflected_total(const uint32_t n, const uint32_t k,
                                const uint32_t d) {
  const uint64_t complement = (uint64_t)k * d - n;
  return complement < n ? (uint32_t)complement : n;
}

uintx compute_bic_no_cache(const uint32_t n, const uint32_t k,
                           const uint32_t d) {
  if (n > INT32_MAX || k > INT32_MAX || d > INT32_MAX) {
    errno = ERANGE;
    return 0;
  }
  if (k == 0 || n == 0 || (uint64_t)n > (uint64_t)k * d) {
    return n == 0 ? 1 : 0;
  }
  const uint32_t total = reflected_total(n, k, d);
  if (!bic_params_supported(total, k, d, BIC_CACHE_NONE)) {
    errno = ERANGE;
    return 0;
  }

  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint32_t j = min(k, total / (d + 1));
  for (uint32_t i = 0; i <= j; ++i) {
    left = compute_bin(k, i, NULL);
    right = compute_bin(total - (d + 1) * i + k - 1, k - 1, NULL);
    compute = (intx)(left * right);
    if (i & 1U) {
      compute = -compute;
    }
    rop += compute;
  }

  return (rop < 0) ? 0 : (uintx)rop;
}

uintx compute_bic_with_sums(const uint32_t n, const uint32_t k,
                            const uint32_t d, intx *partial_sums,
                            const bic_ctx_t ctx) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx compute = 0;
  uintx left = 0;
  uintx right = 0;

  uint32_t j = min(k, n / (d + 1));
  for (uint32_t i = 0; i <= j; ++i) {
    left = ctx->bin(k, i, ctx);
    right = ctx->bin(n - (d + 1) * i + k - 1, k - 1, ctx);
    compute = (intx)(left * right);
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

uintx compute_bic(const uint32_t n, const uint32_t k, const uint32_t d,
                  const bic_ctx_t ctx) {
  if (k == 0 || n == 0 || (uint64_t)n > (uint64_t)k * d) {
    return n == 0 ? 1 : 0;
  }
  return compute_bic_with_sums(reflected_total(n, k, d), k, d, NULL, ctx);
}

uint32_t compute_acc(uintx *rop, const uint32_t n, const uint32_t k,
                     const uint32_t d, const bic_ctx_t ctx) {
  uint32_t i = 0;
  const uint32_t limit = min(n, d);

  rop[0] = 0;
  for (; i <= limit; ++i) {
    rop[i + 1] = ctx->dir(n, k, d, i + 1, ctx);
  }

  return i;
}

uintx compute_dir(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t l, const bic_ctx_t ctx) {
  uint32_t j = min(k, n / (d + 1));
  uint32_t u;

  intx tmp = 0;
  intx rop = 0;

  for (uint32_t i = 0; i <= j; ++i) {
    u = n - (d + 1) * i + k;
    uintx bin_u_l = (u >= l) ? ctx->bin(u - l, k, ctx) : 0;
    tmp = (intx)(ctx->bin(k, i, ctx) * (ctx->bin(u, k, ctx) - bin_u_l));
    if (i & 1U) {
      tmp = -tmp;
    }
    rop += tmp;
  }

  return (rop < 0) ? 0 : (uintx)rop;
}

uintx S_p(const uint32_t n, const uint32_t k, const uint32_t d,
          const uint32_t l, const uint32_t p, const bic_ctx_t ctx) {
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

double exp_total_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                     const uint32_t l, const bic_ctx_t ctx) {
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

double var_total_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                     const uint32_t l, const bic_ctx_t ctx) {
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

double var_single_part(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t l, const bic_ctx_t ctx) {
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

double exp_part_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t j, const uint32_t l, const bic_ctx_t ctx) {
  if (k <= 0) {
    return 0.0;
  }
  return (double)j * exp_total_sum(n, k, d, l, ctx) / k;
}

double var_part_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t j, const uint32_t l, const bic_ctx_t ctx) {
  if (k <= 1) {
    return var_total_sum(n, k, d, l, ctx);
  }

  double lhs =
      ((double)j * (k - j) / (k - 1)) * var_single_part(n, k, d, l, ctx);
  double rhs = ((double)j * (j - 1) / ((double)k * (k - 1))) *
               var_total_sum(n, k, d, l, ctx);

  return dmax(0.0, lhs + rhs);
}

double stddev_part_sum(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t j, const uint32_t l,
                       const bic_ctx_t ctx) {
  return asqrt(var_part_sum(n, k, d, j, l, ctx));
}

void calc_hset(uint32_t k, uint16_t L, uint32_t *H, size_t *h_size) {
  *h_size = 0;
  if (k == 0 || k == 1) {
    return;
  }

  size_t current_size = 0;
  uint32_t last = UINT32_MAX;

  for (size_t p = 1; p <= L; ++p) {
    uint32_t lo = k >> p;
    uint32_t hi = (uint32_t)((k + ((uint32_t)1 << p) - 1) >> p);

    if (hi != last) {
      if (hi > 0) {
        H[current_size++] = hi;
      }
      last = hi;
    }

    if (lo != last) {
      if (lo > 0) {
        H[current_size++] = lo;
      }
      last = lo;
    }
  }

  for (size_t i = 0; i < current_size; ++i) {
    for (size_t j = i + 1; j < current_size; ++j) {
      if (H[i] > H[j]) {
        uint32_t tmp = H[i];
        H[i] = H[j];
        H[j] = tmp;
      }
    }
  }

  *h_size = current_size;
}
