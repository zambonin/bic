#include "utils.h"
#include "math.h"

uint64_t cycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

uint32_t min(const uint32_t a, const uint32_t b) { return (a < b) ? a : b; }

int32_t max(const int32_t a, const int32_t b) { return (a > b) ? a : b; }

void check_valid_bounded_composition(const uint32_t *c, const uint16_t n,
                                     const uint16_t k, const uint16_t d) {
  uint32_t sum = 0;
  for (uint16_t i = 0; i < k; ++i) {
    assert(c[i] <= d);
    sum += c[i];
  }
  assert(sum == n);
}

int cmp(const void *c, const void *d) {
  uintx a = *(uintx *)c;
  uintx b = *(uintx *)d;

  if (a < b) {
    return -1;
  }

  if (a >= b) {
    return 1;
  }

  return 0;
}

static int bsearch_insertion_compare(const void *a, const void *b) {
  bsearch_insertion_state *state = (bsearch_insertion_state *)a;
  state->last_visited = (void *)b;
  return cmp(state->key, b);
}

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width) {
  bsearch_insertion_state state = {key, NULL};
  bsearch(&state, base, nel, width, bsearch_insertion_compare);

  uintx *last = (uintx *)state.last_visited;
  if (*((uintx *)key) < *last) {
    --last;
  }

  return last - (uintx *)base;
}

uint16_t bits_fit_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return 1 + ((uint16_t)lg(alt_inner_bic(n, k, d)));
}

bool bic_geq_2_pow_m(const uint16_t m, const uint16_t n, const uint16_t k,
                     const uint16_t d) {
  return (bool)(alt_inner_bic(n, k, d) >> m);
}

static void linear_search(uint16_t *rop, const uint16_t lo, const uint16_t hi,
                          bool (*p)(const uint16_t, const void *),
                          const void *ctx) {
  uint16_t i;
  for (i = lo; i <= hi && !p(i, ctx); ++i) {
  }
  *rop = i;
}

static bool unimodal(const uint16_t val, const void *ctx) {
  ctx_t *c = (ctx_t *)ctx;
  uint16_t max_n = (c->k * val + 1) / 2;
  return bic_geq_2_pow_m(c->m, max_n, c->k, val);
}

static bool min_n(const uint16_t val, const void *ctx) {
  ctx_t *c = (ctx_t *)ctx;
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

  ctx_t c = {.m = m, .k = k, .v = d};
  linear_search(d, 0, UINT16_MAX, unimodal, &c);
  if (*d == UINT16_MAX) {
    *n = UINT16_MAX;
    return;
  }
  linear_search(n, 0, (k * *d + 1) / 2, min_n, &c);
}

static bool unbounded_parts(const uint16_t val, const void *ctx) {
  ctx_t *c = (ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, val, c->k, val);
}

static bool min_d(const uint16_t val, const void *ctx) {
  ctx_t *c = (ctx_t *)ctx;
  return bic_geq_2_pow_m(c->m, *(c->v), c->k, val);
}

void minver(const uint16_t m, uint16_t *n, const uint16_t k, uint16_t *d) {
  if (k <= 1 || (m >= ((k - 1) * (16 - lg(k - 1))))) {
    return;
  }

  ctx_t c = {.m = m, .k = k, .v = n};
  linear_search(n, 0, UINT16_MAX, unbounded_parts, &c);
  if (*n == UINT16_MAX) {
    *d = UINT16_MAX;
    return;
  }
  linear_search(d, 0, *n, min_d, &c);
}

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d) {
  uint16_t len = bits_fit_bic(n, k, d) / sizeof(uint64_t);
  len += (len == 0);

  uint8_t *message = (uint8_t *)calloc(len, sizeof(uint8_t));
  for (uint16_t i = 0; i < len; ++i) {
    message[i] = random();
  }
  uintx rank = 0;

#if defined(BITINT)
  unsigned char *ptr = (unsigned char *)&rank;
  for (uint16_t i = 0; i < len; ++i) {
    ptr[len - 1 - i] = message[i];
  }
#elif defined(BOOST_FIX_INT) || defined(BOOST_ARB_INT)
  boost::multiprecision::detail::import_bits_fast(rank, message, message + len);
#elif defined(BOOST_MPZ_INT)
  mpz_import(rank.backend().data(), len, 1, sizeof(uint8_t), 0, 0, message);
#elif defined(BOOST_TOM_INT)
  mp_err err =
      mp_unpack(&rank.backend().data(), len, 1, sizeof(uint8_t), 0, 0, message);
  (void)err;
#endif

  free(message);
  return rank % alt_inner_bic(n, k, d);
}
