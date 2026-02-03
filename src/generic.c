#include "generic.h"
#include "api.h"

void bic_generic_unrank(uint32_t *res, uint32_t n, uint32_t k, uint32_t d,
                        const uintx *rank, const bic_iter_strategy_t *strat,
                        const bic_ctx_t ctx) {
  if (k == 1) {
    *res = n;
    return;
  }

  uint32_t k_l, k_r;
  strat->topo(k, &k_l, &k_r);

  bic_iter_state_t state;
  strat->iter_init(&state, n, k_l, k_r, d);

  uint32_t w_l;
  uintx r = *rank;

  while (strat->iter_next(&state, &w_l)) {
    uint32_t w_r = n - w_l;

    uintx vol_l = ctx->comp(w_l, k_l, d, ctx);
    uintx vol_r = ctx->comp(w_r, k_r, d, ctx);
    uintx vol = vol_l * vol_r;

    if (vol == 0)
      continue;

    if (r < vol) {
      uintx r_l, r_r;
      strat->map(r, vol_l, vol_r, w_l, w_r, &r_l, &r_r);

      if (k_l == 1) {
        *res = w_l;
      } else {
        bic_generic_unrank(res, w_l, k_l, d, &r_l, strat, ctx);
      }

      if (k_r == 1) {
        *(res + k_l) = w_r;
      } else {
        bic_generic_unrank(res + k_l, w_r, k_r, d, &r_r, strat, ctx);
      }
      return;
    }
    r -= vol;
  }
}

uintx bic_generic_rank(uint32_t n, uint32_t k, uint32_t d, const uint32_t *res,
                       const bic_iter_strategy_t *strat, const bic_ctx_t ctx) {
  if (k == 1)
    return 0;

  uint32_t k_l, k_r;
  strat->topo(k, &k_l, &k_r);

  uint32_t w_l = 0;
  for (uint32_t i = 0; i < k_l; ++i)
    w_l += res[i];
  uint32_t w_r = n - w_l;

  bic_iter_state_t state;
  strat->iter_init(&state, n, k_l, k_r, d);

  uintx acc = 0;
  uint32_t curr_w;
  while (strat->iter_next(&state, &curr_w)) {
    if (curr_w == w_l) {
      uintx r_l = bic_generic_rank(w_l, k_l, d, res, strat, ctx);
      uintx r_r = bic_generic_rank(w_r, k_r, d, res + k_l, strat, ctx);
      uintx vol_l = ctx->comp(w_l, k_l, d, ctx);
      uintx vol_r = ctx->comp(w_r, k_r, d, ctx);
      return acc + (strat->unmap
                        ? strat->unmap(r_l, r_r, vol_l, vol_r, w_l, w_r)
                        : r_l * vol_r + r_r);
    }
    uint32_t curr_wr = n - curr_w;
    acc += ctx->comp(curr_w, k_l, d, ctx) * ctx->comp(curr_wr, k_r, d, ctx);
  }
  return 0;
}
