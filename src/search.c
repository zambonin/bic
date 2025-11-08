#include "search.h"
#include "math.h"

// https://stackoverflow.com/a/40245287
typedef struct {
  const void *key;
  void *last_visited;
} bsearch_ctx_t;

void linear_search(uint16_t *rop, const uint16_t lo, const uint16_t hi,
                   bool (*p)(const uint16_t, const void *), const void *ctx) {
  uint16_t i;
  for (i = lo; i < hi && !p(i, ctx); ++i) {
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
