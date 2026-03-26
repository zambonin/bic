#include "generic/gray.h"

static void decomp_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void path_gray_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                           uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = n;
  s->val_1 = n - max_wl;
  s->val_2 = n - min_wl;
  s->current_i = s->val_1;
}

static int path_gray_next(bic_path_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_2) {
    return 0;
  }
  *w = s->n - s->current_i++;
  return 1;
}

static void split_gray(uintx rank, uintx vol_l, uintx vol_r, uint32_t w_l,
                       uint32_t w_r, uintx *r_l, uintx *r_r) {
  (void)w_l;
  uintx r_std_l = rank / vol_r;
  uintx r_std_r = rank % vol_r;

  if (w_r & 1U) {
    *r_l = vol_l - 1 - r_std_l;
  } else {
    *r_l = r_std_l;
  }

  *r_r = r_std_r;
}

static uintx unsplit_gray(uintx r_l, uintx r_r, uintx vol_l, uintx vol_r,
                          uint32_t w_l, uint32_t w_r) {
  (void)w_l;
  uintx r_l_used = r_l;
  if (w_r & 1U) {
    r_l_used = vol_l - 1 - r_l;
  }
  return r_l_used * vol_r + r_r;
}

static const bic_path_strategy_t gray_strat = {
    .name = "gray code",
    .decomp = decomp_linear,
    .path_init = path_gray_init,
    .path_next = path_gray_next,
    .split = split_gray,
    .unsplit = unsplit_gray,
};

const bic_path_strategy_t *bic_strat_gray(void) { return &gray_strat; }

void gray_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_gray();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx gray_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_gray();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
