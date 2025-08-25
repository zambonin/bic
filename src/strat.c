#include "strat.h"
#include "math.h"
#include "utils.h"

typedef struct {
  const uint16_t m;
  const uint16_t k;
  const uint16_t *v;
} param_ctx_t;

bool bic_geq_2_pow_m(const uint16_t m, const uint16_t n, const uint16_t k,
                     const uint16_t d) {
  return (bool)(alt_compute_bic(n, k, d) >> m);
}

bool unimodal(const uint16_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  uint16_t max_n = (c->k * val + 1) / 2;
  return bic_geq_2_pow_m(c->m, max_n, c->k, val);
}

bool min_n(const uint16_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, val, c->k, *(c->v));
}

void mingen(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d) {
  /*
   * |C(n, k, d)| is bounded by |C(n, k, n)| = \binom{n + k - 1}{k - 1}.
   *
   * The binomial is approximately n^{k - 1} / (k - 1)!. Taking the base-2 log,
   * we have (k - 1) * (lg(n) - lg(k - 1)), ignoring the smaller part of the
   * denominator, and we can set lg(n) = 16 for uint16_t.
   */
  if (k <= 1 || (m >= ((k - 1) * (16 - lg(k - 1))))) {
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = d};
  linear_search(d, 0, UINT16_MAX, unimodal, &c);
  if (*d == UINT16_MAX) {
    *n = UINT16_MAX;
    return;
  }
  linear_search(n, 0, (k * *d + 1) / 2, min_n, &c);
}

bool unbounded_parts(const uint16_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, val, c->k, val);
}

bool min_d(const uint16_t val, const void *ctx) {
  param_ctx_t *c = (param_ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, *(c->v), c->k, val);
}

void minver(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d) {
  if (k <= 1 || (m >= ((k - 1) * (16 - lg(k - 1))))) {
    return;
  }

  param_ctx_t c = {.m = m, .k = k, .v = n};
  linear_search(n, 0, UINT16_MAX, unbounded_parts, &c);
  if (*n == UINT16_MAX) {
    *d = UINT16_MAX;
    return;
  }
  linear_search(d, 0, *n, min_d, &c);
}

void gen_params_random(const uint16_t m, uint16_t *n, const uint16_t k,
                       uint16_t *d) {
  do {
    *d = (random() % 64) + 1;
    *n = (random() % ((k * *d + 1) / 2)) + 1;
  } while (!bic_geq_2_pow_m(*n, k, *d, m));
}
