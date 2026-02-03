#include "spiral.h"
#include "cache.h"
#include "common.h"
#include "generic.h"
#include "math.h"

void spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                   const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t p = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; rop[i] = p, --i, it_n -= p) {
    uint32_t min_p = 0;
    uintx max_rem = (uintx)i * d;
    if ((uintx)it_n > max_rem) {
      min_p = it_n - (uint32_t)max_rem;
    }
    uint32_t max_p = (it_n < d) ? it_n : d;

    uint32_t start_p = it_n / (i + 1);
    if (start_p < min_p)
      start_p = min_p;
    if (start_p > max_p)
      start_p = max_p;

    uint32_t range_up = max_p - start_p;
    uint32_t range_down = start_p - min_p;
    uint32_t max_range = (range_up > range_down) ? range_up : range_down;

    for (uint32_t step = 0;; ++step) {
      if (step == 0) {
        p = start_p;
      } else {
        uint16_t mag = (step + 1) / 2;
        if (mag > max_range + 1)
          break;
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
    uintx max_rem = (uintx)i * d;
    if ((uintx)it_n > max_rem) {
      min_p = it_n - (uint32_t)max_rem;
    }
    uint32_t max_p = (it_n < d) ? it_n : d;

    uint32_t start_p = it_n / (i + 1);
    if (start_p < min_p)
      start_p = min_p;
    if (start_p > max_p)
      start_p = max_p;

    uint32_t range_up = max_p - start_p;
    uint32_t range_down = start_p - min_p;
    uint32_t max_range = (range_up > range_down) ? range_up : range_down;

    for (uint32_t step = 0;; ++step) {
      if (step == 0) {
        p = start_p;
      } else {
        uint16_t mag = (step + 1) / 2;
        if (mag > max_range + 1)
          break;
        if (step % 2 != 0) {
          p = start_p + mag;
        } else {
          p = start_p - mag;
        }
      }

      if (p >= min_p && p <= max_p) {
        if (p == target)
          break;
        count = ctx->comp(it_n - p, i, d, ctx);
        rank += count;
      }
    }
  }

  return rank;
}

static void topo_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void iter_spiral_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                             uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  uint32_t min_wr = n - max_wl;
  uint32_t max_wr = n - min_wl;

  uint32_t mean_wr = (uint32_t)((uintx)n * k_r / (k_l + k_r));
  if (mean_wr < min_wr)
    mean_wr = min_wr;
  if (mean_wr > max_wr)
    mean_wr = max_wr;

  s->n = n;
  s->val_1 = min_wr;
  s->val_2 = max_wr;
  s->k_l = mean_wr;
  s->current_i = 0;
}

static int iter_spiral_next(bic_iter_state_t *s, uint32_t *w) {
  uint32_t center = s->k_l;
  uint32_t min_wr = s->val_1;
  uint32_t max_wr = s->val_2;

  while (1) {
    uint32_t step = s->current_i++;
    int16_t offset =
        (step == 0) ? 0 : ((step & 1) ? (step + 1) / 2 : -(step / 2));
    int32_t wr = (int32_t)center + offset;

    if (offset > 0 && (int64_t)wr > (int64_t)max_wr &&
        (int64_t)center - offset < (int64_t)min_wr)
      return 0;
    if (offset < 0 && (int64_t)wr < (int64_t)min_wr &&
        (int64_t)center - offset > (int64_t)max_wr)
      return 0;

    if ((int64_t)wr >= (int64_t)min_wr && (int64_t)wr <= (int64_t)max_wr) {
      *w = s->n - (uint32_t)wr;
      return 1;
    }
  }
}

static const bic_iter_strategy_t spiral_strat = {
    .name = "Spiral",
    .topo = topo_linear,
    .iter_init = iter_spiral_init,
    .iter_next = iter_spiral_next,
    .map = map_std,
    .unmap = unmap_std,
};

const bic_iter_strategy_t *bic_strat_spiral(void) { return &spiral_strat; }

void spiral_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                           const uint32_t d, const uintx r,
                           const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_spiral();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx spiral_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                          const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_spiral();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
