#include "colex.h"
#include "cache.h"
#include "common.h"
#include "math.h"
#include "search.h"

void colex_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                  const uint16_t d, const uintx r, const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0; count = ctx->comp(it_n - part, i, d, ctx), rank >= count;
         ++part, rank -= count) {
    }
  }

  rop[0] = it_n;
}

void colex_unrank_part_sums(uint32_t *rop, const uint16_t n, const uint16_t k,
                            const uint16_t d, const uintx r,
                            const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;

  uint16_t j = min(k - 1, it_n / (d + 1));
  intx *prev_sum = intx_alloc(j + 1);

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    intx left = 0;
    intx right = compute_bic_with_sums(it_n, i, d, prev_sum, ctx);

    for (part = 0; rank >= (uintx)right; ++part) {
      left = right;
      uint16_t numerator = (it_n - part);
      uint16_t denominator = (it_n - part) + i - 1;

      j = min(i, (it_n - part - 1) / (d + 1));
      for (uint16_t m = 0; m <= j; ++m) {
        prev_sum[m] *= numerator;
        prev_sum[m] /= denominator;
        right += prev_sum[m];
        numerator -= (d + 1);
        denominator -= (d + 1);
      }
    }

    rank -= (uintx)left;
  }

  rop[0] = it_n;

  intx_free(prev_sum);
}

void colex_unrank_acc_linear(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;
  uintx count = 0;
  uintx *sums = uintx_alloc(d + 3);

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    (void)ctx->acc(sums, it_n, i, d, ctx);
    for (part = 0; count = sums[part + 1], rank >= count; ++part) {
    }
    rank -= sums[part];
  }

  rop[0] = it_n;
  uintx_free(sums);
}

void colex_unrank_acc_bisect(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;
  uintx *sums = uintx_alloc(d + 3);

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    size_t length = ctx->acc(sums, it_n, i, d, ctx);
    part = bsearch_insertion(&rank, sums, length, sizeof(uintx));
    rank -= sums[part];
  }

  rop[0] = it_n;
  uintx_free(sums);
}

void colex_unrank_acc_direct(uint32_t *rop, const uint16_t n, const uint16_t k,
                             const uint16_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    part = 0;
    for (uint16_t c = min(it_n, d); c > 0;) {
      uint16_t step = (c / 2) + 1;
      count = ctx->dir(it_n, i, d, part + step, ctx);
      if (rank >= count) {
        part += step;
        c -= step;
      } else {
        c = step - 1;
      }
    }
    rank -= ctx->dir(it_n, i, d, part, ctx);
  }

  rop[0] = it_n;
}

uintx colex_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                 const uint32_t *comb, const bic_ctx_t ctx) {
  uintx rank = 0;
  uint16_t it_n = n;

  for (uint16_t i = k - 1; i > 0; it_n -= comb[i], --i) {
    for (uint16_t j = 0; j < comb[i];
         rank += ctx->comp(it_n - j, i, d, ctx), ++j) {
    }
  }

  return rank;
}
