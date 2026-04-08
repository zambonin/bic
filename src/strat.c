#include "strat.h"
#include "math.h"
#include "search.h"

typedef struct {
  const uint32_t m;
  const uint32_t k;
  const uint32_t *v;
} param_ctx_t;

bool bic_geq_2_pow_m(const uint32_t m, const uint32_t n, const uint32_t k,
                     const uint32_t d) {
#if defined(BITINT) || defined(BOOST_FIX_INT)
  if (m >= BIT_LENGTH - 1) {
    return false;
  }
#endif
  return (bool)(compute_bic_no_cache(n, k, d) >> m);
}

bool unimodal(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  uint64_t max_n_long = ((uint64_t)c->k * val + 1) / 2;
  uint32_t max_n =
      (max_n_long > UINT32_MAX) ? UINT32_MAX : (uint32_t)max_n_long;
  return bic_geq_2_pow_m(c->m, max_n, c->k, val);
}

bool min_n(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, val, c->k, *(c->v));
}

void mingen(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d) {
  if (k <= 1 || (m >= ((k - 1) * (32 - lg(k - 1))))) {
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = d};
  exponential_search(d, 0, UINT32_MAX, unimodal, &c);
  if (*d == UINT32_MAX) {
    *n = UINT32_MAX;
    return;
  }
  uint64_t hi_long = ((uint64_t)k * *d + 1) / 2;
  exponential_search(n, 0, (uint32_t)hi_long, min_n, &c);
}

bool unbounded_parts(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, val, c->k, val);
}

bool min_d(const uint32_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, *(c->v), c->k, val);
}

void minver(const uint32_t m, uint32_t *n, const uint32_t k, uint32_t *d) {
  if (k <= 1 || (m >= ((k - 1) * (32 - lg(k - 1))))) {
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = n};
  exponential_search(n, 0, UINT32_MAX, unbounded_parts, &c);
  if (*n == UINT32_MAX) {
    *d = UINT32_MAX;
    return;
  }
  exponential_search(d, 0, *n, min_d, &c);
}

void gen_params_random(const uint32_t m, uint32_t *n, const uint32_t k,
                       uint32_t *d) {
  uint32_t chosen_d = 0;
  uint32_t choices = 0;

  for (uint32_t candidate = 1; candidate <= 64; ++candidate) {
    const uint64_t center = ((uint64_t)k * candidate + 1) / 2;
    if (center == 0 || center > INT32_MAX ||
        !bic_params_supported((uint32_t)center, k, candidate, BIC_CACHE_NONE) ||
        !bic_geq_2_pow_m(m, (uint32_t)center, k, candidate)) {
      continue;
    }
    ++choices;
    if (random() % choices == 0) {
      chosen_d = candidate;
    }
  }

  if (choices == 0) {
    *n = *d = 0;
    return;
  }

  const uint32_t center = (uint32_t)(((uint64_t)k * chosen_d + 1) / 2);
  uint32_t lo = 1;
  uint32_t hi = center;
  while (lo < hi) {
    const uint32_t mid = lo + (hi - lo) / 2;
    if (bic_geq_2_pow_m(m, mid, k, chosen_d)) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }

  *d = chosen_d;
  *n = lo + (uint32_t)(random() % (center - lo + 1));
}
