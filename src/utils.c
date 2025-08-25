#include "utils.h"
#include "math.h"

// https://stackoverflow.com/a/40245287
typedef struct {
  const void *key;
  void *last_visited;
} bsearch_ctx_t;

void linear_search(uint16_t *rop, const uint16_t lo, const uint16_t hi,
                   bool (*p)(const uint16_t, const void *), const void *ctx) {
  uint16_t i;
  for (i = lo; i <= hi && !p(i, ctx); ++i) {
  }
  *rop = i;
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

int bsearch_insertion_compare(const void *a, const void *b) {
  bsearch_ctx_t *state = (bsearch_ctx_t *)a;
  state->last_visited = (void *)b;
  return cmp(state->key, b);
}

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width) {
  bsearch_ctx_t state = {key, NULL};
  bsearch(&state, base, nel, width, bsearch_insertion_compare);

  uintx *last = (uintx *)state.last_visited;
  if (*((uintx *)key) < *last) {
    --last;
  }

  return last - (uintx *)base;
}

uint16_t bits_fit_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return 1 + ((uint16_t)lg(alt_compute_bic(n, k, d)));
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
  return rank % alt_compute_bic(n, k, d);
}
