#include "search.h"
#include "math.h"

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

size_t uintx_upper_bound(const uintx *values, size_t count, const uintx *key) {
  size_t lo = 0, hi = count;
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (values[mid] <= *key) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}
