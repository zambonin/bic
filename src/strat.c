#include "strat.h"
#include "math.h"
#include "search.h"
#include <errno.h>

typedef struct {
  const uint32_t m;
  const uint32_t k;
  const uint32_t *v;
} param_ctx_t;

bic_count_comparison_t bic_compare_2_pow_m(const uint32_t m, const uint32_t n,
                                           const uint32_t k, const uint32_t d) {
  uintx count;
  if (!bic_count_checked(&count, n, k, d)) {
    return BIC_COUNT_UNSUPPORTED;
  }
#if defined(BITINT) || defined(BOOST_FIX_INT)

  if (m >= BIT_LENGTH - 1) {
    return BIC_COUNT_BELOW;
  }
#endif
  return (count >> m) != 0 ? BIC_COUNT_AT_LEAST : BIC_COUNT_BELOW;
}

bool bic_geq_2_pow_m(const uint32_t m, const uint32_t n, const uint32_t k,
                     const uint32_t d) {
  return bic_compare_2_pow_m(m, n, k, d) == BIC_COUNT_AT_LEAST;
}

static void search_failed(uint32_t *n, uint32_t *d) {
  *n = *d = UINT32_MAX;
  errno = ERANGE;
}

static bool eligible(const uint32_t m, const uint32_t n, const uint32_t k,
                     const uint32_t d) {
  if (!bic_params_supported(n, k, d, BIC_CACHE_NONE)) {
    return false;
  }
  return bic_compare_2_pow_m(m, n, k, d) == BIC_COUNT_AT_LEAST;
}

static uint32_t largest_supported_n(const uint32_t upper, const uint32_t k,
                                    const uint32_t d, const bool unbounded) {

  uint32_t lo = 0;
  uint32_t hi = upper;
  while (lo < hi) {
    const uint32_t mid = lo + (hi - lo) / 2 + (hi - lo) % 2;
    if (bic_params_supported(mid, k, unbounded ? mid : d, BIC_CACHE_NONE)) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  return lo;
}

bool unimodal(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  if (val > INT32_MAX || c->k > INT32_MAX) {
    return false;
  }
  uint64_t max_n_long = ((uint64_t)c->k * val + 1) / 2;
  uint32_t max_n = (max_n_long > INT32_MAX) ? INT32_MAX : (uint32_t)max_n_long;
  max_n = largest_supported_n(max_n, c->k, val, false);
  return eligible(c->m, max_n, c->k, val);
}

bool min_n(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return eligible(c->m, val, c->k, *(c->v));
}

void mingen(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d) {
  if (m == 0 && k <= INT32_MAX) {
    *n = *d = 0;
    errno = 0;
    return;
  }
  if (k <= 1 || k > INT32_MAX) {
    search_failed(n, d);
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = d};
  exponential_search(d, 0, INT32_MAX, unimodal, &c);
  if (!unimodal(*d, &c)) {
    search_failed(n, d);
    return;
  }
  uint64_t hi_long = ((uint64_t)k * *d + 1) / 2;
  const uint32_t hi = largest_supported_n(
      hi_long > INT32_MAX ? INT32_MAX : (uint32_t)hi_long, k, *d, false);
  binary_search_first_true(n, 0, hi, min_n, &c);
  errno = 0;
}

bool unbounded_parts(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return eligible(c->m, val, c->k, val);
}

bool min_d(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return eligible(c->m, *(c->v), c->k, val);
}

void minver(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d) {
  if (m == 0 && k <= INT32_MAX) {
    *n = *d = 0;
    errno = 0;
    return;
  }
  if (k <= 1 || k > INT32_MAX) {
    search_failed(n, d);
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = n};
  const uint32_t hi = largest_supported_n(INT32_MAX, k, 0, true);
  if (!eligible(m, hi, k, hi)) {
    search_failed(n, d);
    return;
  }
  binary_search_first_true(n, 0, hi, unbounded_parts, &c);
  binary_search_first_true(d, 0, *n, min_d, &c);
  errno = 0;
}

void gen_params_random(const uint32_t m, uint32_t *n, const uint32_t k,
                       uint32_t *d) {
  if (m == 0 && k <= INT32_MAX) {
    *n = *d = 0;
    errno = 0;
    return;
  }
  uint32_t chosen_d = 0;
  uint32_t choices = 0;

  for (uint32_t candidate = 1; candidate <= 64; ++candidate) {
    const uint64_t center = ((uint64_t)k * candidate + 1) / 2;
    if (center == 0 || center > INT32_MAX ||
        !bic_params_supported((uint32_t)center, k, candidate, BIC_CACHE_NONE) ||
        !eligible(m, (uint32_t)center, k, candidate)) {
      continue;
    }
    ++choices;
    if (random() % choices == 0) {
      chosen_d = candidate;
    }
  }

  if (choices == 0) {
    search_failed(n, d);
    return;
  }

  const uint32_t center = (uint32_t)(((uint64_t)k * chosen_d + 1) / 2);
  uint32_t lo = 1;
  uint32_t hi = center;
  while (lo < hi) {
    const uint32_t mid = lo + (hi - lo) / 2;
    if (eligible(m, mid, k, chosen_d)) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }

  *d = chosen_d;
  *n = lo + (uint32_t)(random() % (center - lo + 1));
  errno = 0;
}
