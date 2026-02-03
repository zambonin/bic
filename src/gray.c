#include "gray.h"
#include "generic.h"
#include "math.h"
#include "types.h"

void gray_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                 const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = r;
  uint32_t part = 0;
  uintx count = 0;

  for (uint32_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0; count = ctx->comp(it_n - part, i, d, ctx), rank >= count;
         ++part, rank -= count) {
    }
    if (part & 1U) {
      rank = count - 1 - rank;
    }
  }

  rop[0] = it_n;
}

uintx gray_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                const uint32_t *comb, const bic_ctx_t ctx) {
  uint32_t it_n = n;
  uintx rank = 0;
  uint32_t part = 0;

  for (uint32_t i = k - 1, p = n & 1U, parity = 0;
       part = comb[i], parity = it_n & 1U, i > 0; --i, it_n -= part) {
    if (parity == p) {
      for (uint32_t j = 0; j < part; ++j) {
        rank += ctx->comp(it_n - j, i, d, ctx);
      }
    } else {
      for (uint32_t j = part + 1; j < min(it_n, d) + 1; ++j) {
        rank += ctx->comp(it_n - j, i, d, ctx);
      }
    }
  }

  return rank;
}

static void topo_linear(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k - 1;
  *k_r = 1;
}

static void iter_gray_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                           uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = n;
  s->val_1 = n - max_wl;
  s->val_2 = n - min_wl;
  s->current_i = s->val_1;
}

static int iter_gray_next(bic_iter_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_2)
    return 0;
  *w = s->n - s->current_i++;
  return 1;
}

static void map_gray(uintx rank, uintx vol_l, uintx vol_r, uint32_t w_l,
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

static uintx unmap_gray(uintx r_l, uintx r_r, uintx vol_l, uintx vol_r,
                        uint32_t w_l, uint32_t w_r) {
  (void)w_l;
  uintx r_l_used = r_l;
  if (w_r & 1U) {
    r_l_used = vol_l - 1 - r_l;
  }
  return r_l_used * vol_r + r_r;
}

static const bic_iter_strategy_t gray_strat = {
    .name = "Gray",
    .topo = topo_linear,
    .iter_init = iter_gray_init,
    .iter_next = iter_gray_next,
    .map = map_gray,
    .unmap = unmap_gray,
};

const bic_iter_strategy_t *bic_strat_gray(void) { return &gray_strat; }

void gray_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_gray();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx gray_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_gray();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
