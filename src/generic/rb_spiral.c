#include "generic/rb_spiral.h"

static void decomp_balanced(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k / 2;
  *k_r = k - *k_l;
}

static void path_rb_spiral_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                                uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = min_wl;
  s->val_1 = max_wl;
  s->val_2 = (min_wl + max_wl) / 2;
  s->current_i = 0;
}

static int path_rb_spiral_next(bic_path_state_t *s, uint32_t *w) {
  uint32_t min_wl = s->n;
  uint32_t max_wl = s->val_1;
  uint32_t mid = s->val_2;

  while (1) {
    uint32_t jump = (s->current_i + 1) / 2;
    int sign = s->current_i % 2;
    uint32_t val;

    if (sign) {
      if (jump > mid) {
        val = min_wl - 1;
      } else {
        val = mid - jump;
      }
    } else {
      val = mid + jump;
    }

    s->current_i++;

    if (val >= min_wl && val <= max_wl) {
      *w = val;
      return 1;
    }

    if (jump > (max_wl > mid ? max_wl - mid : 0) &&
        jump > (mid > min_wl ? mid - min_wl : 0)) {
      return 0;
    }
  }
}

static const bic_path_strategy_t rb_spiral_strat = {
    .name = "rb-spiral",
    .decomp = decomp_balanced,
    .path_init = path_rb_spiral_init,
    .path_next = path_rb_spiral_next,
    .split = split_std,
    .unsplit = unsplit_std,
};

const bic_path_strategy_t *bic_strat_rb_spiral(void) {
  return &rb_spiral_strat;
}

void rb_spiral_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                              const uint32_t d, const uintx r,
                              const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_rb_spiral();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx rb_spiral_generic_rank(const uint32_t n, const uint32_t k,
                             const uint32_t d, const uint32_t *comb,
                             const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_rb_spiral();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
