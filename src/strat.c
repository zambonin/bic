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
  do {
    *d = (random() % 64) + 1;
    *n = (random() % ((k * *d + 1) / 2)) + 1;
  } while (!bic_geq_2_pow_m(m, *n, k, *d));
}

void map_std(uintx rank, uintx vol_l, uintx vol_r, uint32_t w_l, uint32_t w_r,
             uintx *r_l, uintx *r_r) {
  (void)w_l;
  (void)w_r;
  (void)vol_l;
  *r_l = rank / vol_r;
  *r_r = rank % vol_r;
}

uintx unmap_std(uintx r_l, uintx r_r, uintx vol_l, uintx vol_r, uint32_t w_l,
                uint32_t w_r) {
  (void)w_l;
  (void)w_r;
  (void)vol_l;
  return r_l * vol_r + r_r;
}
