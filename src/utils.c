#include "utils.h"

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
