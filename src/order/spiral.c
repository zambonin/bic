#include "order/spiral.h"
#include "common.h"
#include "math.h"
#include "types.h"

void spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                   const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t p = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; rop[i] = p, --i, it_n -= p) {
    uint32_t min_p = 0;
    uint64_t max_rem = (uint64_t)i * d;
    if ((uint64_t)it_n > max_rem) {
      min_p = it_n - (uint32_t)max_rem;
    }
    uint32_t max_p = (it_n < d) ? it_n : d;

    uint32_t start_p = it_n / (i + 1);
    if (start_p < min_p) {
      start_p = min_p;
    }
    if (start_p > max_p) {
      start_p = max_p;
    }

    uint32_t range_up = max_p - start_p;
    uint32_t range_down = start_p - min_p;
    uint32_t max_range = (range_up > range_down) ? range_up : range_down;

    for (uint32_t step = 0;; ++step) {
      if (step == 0) {
        p = start_p;
      } else {
        uint32_t mag = (step + 1) / 2;
        if (mag > max_range + 1) {
          break;
        }
        if (step % 2 != 0) {
          p = start_p + mag;
        } else {
          p = start_p - mag;
        }
      }

      if (p >= min_p && p <= max_p) {
        count = ctx->comp(it_n - p, i, d, ctx);
        if (rank < count) {
          break;
        }
        rank -= count;
      }
    }
  }

  rop[0] = it_n;
}

uintx spiral_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t *comb, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = 0;
  uint32_t p = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; it_n -= comb[i], --i) {
    uint32_t target = comb[i];
    uint32_t min_p = 0;
    uint64_t max_rem = (uint64_t)i * d;
    if ((uint64_t)it_n > max_rem) {
      min_p = it_n - (uint32_t)max_rem;
    }
    uint32_t max_p = (it_n < d) ? it_n : d;

    uint32_t start_p = it_n / (i + 1);
    if (start_p < min_p) {
      start_p = min_p;
    }
    if (start_p > max_p) {
      start_p = max_p;
    }

    uint32_t range_up = max_p - start_p;
    uint32_t range_down = start_p - min_p;
    uint32_t max_range = (range_up > range_down) ? range_up : range_down;

    for (uint32_t step = 0;; ++step) {
      if (step == 0) {
        p = start_p;
      } else {
        uint32_t mag = (step + 1) / 2;
        if (mag > max_range + 1) {
          break;
        }
        if (step % 2 != 0) {
          p = start_p + mag;
        } else {
          p = start_p - mag;
        }
      }

      if (p >= min_p && p <= max_p) {
        if (p == target) {
          break;
        }
        count = ctx->comp(it_n - p, i, d, ctx);
        rank += count;
      }
    }
  }

  return rank;
}
