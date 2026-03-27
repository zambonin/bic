#include "generic/rbo.h"

static void decomp_balanced(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k / 2;
  *k_r = k - *k_l;
}

static void path_rbo_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                          uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->current_i = min_wl;
  s->val_1 = max_wl;
}

static int path_rbo_next(bic_path_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_1) {
    return 0;
  }
  *w = s->current_i++;
  return 1;
}

static const bic_path_strategy_t rbo_strat = {
    .name = "RBO",
    .decomp = decomp_balanced,
    .path_init = path_rbo_init,
    .path_next = path_rbo_next,
    .split = split_std,
    .unsplit = unsplit_std,
};

const bic_path_strategy_t *bic_strat_rbo(void) { return &rbo_strat; }

void rbo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                        const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_rbo();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx rbo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_rbo();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
