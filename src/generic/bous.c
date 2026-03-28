#include "generic/bous.h"

static void decomp_balanced(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k / 2;
  *k_r = k - *k_l;
}

static void path_bous_init(bic_path_state_t *s, uint32_t n, uint32_t k_l,
                           uint32_t k_r, uint32_t d) {
  (void)k_l;
  (void)k_r;
  (void)d;

  s->n = n;
  s->val_1 = 0;
  s->val_2 = n / 2;
  s->current_i = 0;
}

static int path_bous_next(bic_path_state_t *s, uint32_t *w) {
  uint32_t j = s->val_1;
  uint32_t max_j = s->val_2;
  uint32_t n = s->n;
  uint32_t phase = s->current_i;

  if (j > max_j) {
    return 0;
  }

  uint32_t candidate;

  if (phase == 0) {
    candidate = j;
    if (j == n - j) {
      s->val_1++;
      s->current_i = 0;
    } else {
      s->current_i = 1;
    }
  } else {
    candidate = n - j;
    s->val_1++;
    s->current_i = 0;
  }

  *w = candidate;
  return 1;
}

static const bic_path_strategy_t bous_strat = {
    .name = "Bous",
    .decomp = decomp_balanced,
    .path_init = path_bous_init,
    .path_next = path_bous_next,
    .split = split_std,
    .unsplit = unsplit_std,
};

const bic_path_strategy_t *bic_strat_bous(void) { return &bous_strat; }

void bous_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_bous();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx bous_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_path_strategy_t *strat = bic_strat_bous();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
