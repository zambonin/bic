#include "generic/colex.h"

static void decomp_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void path_colex_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                            uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = n;
  s->val_1 = n - max_wl;
  s->val_2 = n - min_wl;
  s->current_i = s->val_1;
}

static int path_colex_next(bic_path_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_2) {
    return 0;
  }
  *w = s->n - s->current_i++;
  return 1;
}

static const bic_path_strategy_t colex_strat = {
    .name = "colexicographic",
    .decomp = decomp_linear,
    .path_init = path_colex_init,
    .path_next = path_colex_next,
    .split = split_std,
    .unsplit = unsplit_std,
};

const bic_path_strategy_t *bic_strat_colex(void) { return &colex_strat; }

void colex_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                          const uint32_t d, const uintx r,
                          const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_colex();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx colex_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                         const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_colex();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
