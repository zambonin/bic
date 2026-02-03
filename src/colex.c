#include "colex.h"
#include "cache.h"
#include "common.h"
#include "generic.h"
#include "math.h"
#include "search.h"

void colex_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                  const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0; count = ctx->comp(it_n - part, i, d, ctx), rank >= count;
         ++part, rank -= count) {
    }
  }

  rop[0] = it_n;
}

void colex_unrank_part_sums(uint32_t *rop, const uint32_t n, const uint32_t k,
                            const uint32_t d, const uintx r,
                            const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;

  uint32_t j = min(k - 1, it_n / (d + 1));
  intx *prev_sum = intx_alloc(j + 1);

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    intx left = 0;
    intx right = compute_bic_with_sums(it_n, i, d, prev_sum, ctx);

    for (part = 0; rank >= (uintx)right; ++part) {
      left = right;
      uint32_t numerator = (it_n - part);
      uint32_t denominator = (it_n - part) + i - 1;

      j = min(i, (it_n - part - 1) / (d + 1));
      for (uint32_t m = 0; m <= j; ++m) {
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

void colex_unrank_acc_linear(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;
  uintx count = 0;
  uintx *sums = uintx_alloc(d + 2);

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    (void)ctx->acc(sums, it_n, i, d, ctx);
    for (part = 0; count = sums[part + 1], rank >= count; ++part) {
    }
    rank -= sums[part];
  }

  rop[0] = it_n;
  uintx_free(sums);
}

void colex_unrank_acc_bisect(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;
  uintx *sums = uintx_alloc(d + 2);

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    size_t length = ctx->acc(sums, it_n, i, d, ctx);
    part = bsearch_insertion(&rank, sums, length, sizeof(uintx));
    rank -= sums[part];
  }

  rop[0] = it_n;
  uintx_free(sums);
}

void colex_unrank_acc_direct(uint32_t *rop, const uint32_t n, const uint32_t k,
                             const uint32_t d, const uintx r,
                             const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    part = 0;
    for (uint32_t c = min(it_n, d); c > 0;) {
      uint32_t step = (c / 2) + 1;
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

uintx colex_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                 const uint32_t *comb, const bic_ctx_t ctx) {
  uintx rank = 0;
  uint32_t it_n = n;

  for (uint32_t i = k - 1; i > 0; it_n -= comb[i], --i) {
    for (uint32_t j = 0; j < comb[i];
         rank += ctx->comp(it_n - j, i, d, ctx), ++j) {
    }
  }

  return rank;
}

uintx colex_rank_acc_direct(const uint32_t n, const uint32_t k,
                            const uint32_t d, const uint32_t *comb,
                            const bic_ctx_t ctx) {
  uintx rank = 0;
  uint32_t it_n = n;

  for (uint32_t i = k - 1; i > 0; it_n -= comb[i], --i) {
    rank += ctx->dir(it_n, i, d, comb[i], ctx);
  }

  return rank;
}

static void topo_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void iter_colex_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                            uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = n;
  s->val_1 = n - max_wl;
  s->val_2 = n - min_wl;
  s->current_i = s->val_1;
}

static int iter_colex_next(bic_iter_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_2)
    return 0;
  *w = s->n - s->current_i++;
  return 1;
}

static const bic_iter_strategy_t colex_strat = {
    .name = "Colex",
    .topo = topo_linear,
    .iter_init = iter_colex_init,
    .iter_next = iter_colex_next,
    .map = map_std,
    .unmap = unmap_std,
};

const bic_iter_strategy_t *bic_strat_colex(void) { return &colex_strat; }

void colex_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                          const uint32_t d, const uintx r,
                          const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_colex();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx colex_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                         const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_colex();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
