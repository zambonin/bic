#include "search.h"
#include "math.h"
#include <stdlib.h>

typedef struct {
  const void *key;
  void *last_visited;
} bsearch_ctx_t;

void linear_search(uint32_t *rop, const uint32_t lo, const uint32_t hi,
                   bool (*p)(const uint32_t, const void *), const void *ctx) {
  uint32_t i;
  for (i = lo; i < hi && !p(i, ctx); ++i) {
  }
  *rop = i;
}

void binary_search_first_true(uint32_t *rop, const uint32_t lo,
                              const uint32_t hi,
                              bool (*p)(const uint32_t, const void *),
                              const void *ctx) {
  uint32_t l = lo;
  uint32_t h = hi;

  while (l < h) {
    uint32_t mid = l + (h - l) / 2;
    if (p(mid, ctx)) {
      h = mid;
    } else {
      l = mid + 1;
    }
  }
  *rop = l;
}

void exponential_search(uint32_t *rop, const uint32_t lo, const uint32_t hi,
                        bool (*p)(const uint32_t, const void *),
                        const void *ctx) {
  if (p(lo, ctx)) {
    *rop = lo;
    return;
  }

  uint32_t bound = 1;
  while (bound < hi - lo && !p(lo + bound, ctx)) {
    if (bound > (hi - lo) / 2) {
      bound = hi - lo;
      break;
    }
    bound *= 2;
  }

  uint32_t start = lo + bound / 2;
  uint32_t end = (lo + bound < hi) ? (lo + bound) : hi;
  binary_search_first_true(rop, start, end, p, ctx);
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
