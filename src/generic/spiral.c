#include "generic/spiral.h"

static void decomp_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void path_spiral_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                             uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  uint32_t min_wr = n - max_wl;
  uint32_t max_wr = n - min_wl;

  uint32_t mean_wr = (uint32_t)((uint64_t)n * k_r / (k_l + k_r));
  if (mean_wr < min_wr) {
    mean_wr = min_wr;
  }
  if (mean_wr > max_wr) {
    mean_wr = max_wr;
  }

  s->n = n;
  s->val_1 = min_wr;
  s->val_2 = max_wr;
  s->k_l = mean_wr;
  s->current_i = 0;
}

static int path_spiral_next(bic_path_state_t *s, uint32_t *w) {
  uint32_t center = s->k_l;
  uint32_t min_wr = s->val_1;
  uint32_t max_wr = s->val_2;

  while (1) {
    uint32_t step = s->current_i++;
    int64_t offset =
        (step == 0) ? 0 : ((step & 1) ? (step + 1) / 2 : -(int64_t)(step / 2));
    int64_t wr = (int64_t)center + offset;

    if (offset > 0 && wr > max_wr && (int64_t)center - offset < min_wr) {
      return 0;
    }
    if (offset < 0 && wr < min_wr && (int64_t)center - offset > max_wr) {
      return 0;
    }

    if (wr >= min_wr && wr <= max_wr) {
      *w = s->n - (uint32_t)wr;
      return 1;
    }
  }
}

static const bic_path_strategy_t spiral_strat = {
    .name = "spiral",
    .decomp = decomp_linear,
    .path_init = path_spiral_init,
    .path_next = path_spiral_next,
    .split = split_std,
    .unsplit = unsplit_std,
};

const bic_path_strategy_t *bic_strat_spiral(void) { return &spiral_strat; }

void spiral_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                           const uint32_t d, const uintx r,
                           const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_spiral();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx spiral_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                          const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_spiral();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
