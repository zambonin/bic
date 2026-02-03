#include "rbo.h"
#include "generic.h"
#include "math.h"
#include "types.h"

void inner_rbo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                      const uint32_t d, const uintx r, const uint32_t start,
                      const bic_ctx_t ctx) {
  uintx rank = r;

  if (k == 1) {
    rop[start] = n;
    return;
  }

  uint32_t left = (uint32_t)(k / 2);
  uint32_t right = k - left;

  uint32_t leftSum = 0;
  uint32_t rightSum = 0;
  uintx rightPoints = 0;

  for (uintx count = 0; leftSum <= min(n, left * d); ++leftSum, rank -= count) {
    rightSum = n - leftSum;
    rightPoints = ctx->comp(rightSum, right, d, ctx);
    count = ctx->comp(leftSum, left, d, ctx) * rightPoints;
    if (rank < count) {
      break;
    }
  }

  uintx leftRank = rank / rightPoints;
  uintx rightRank = rank % rightPoints;

  inner_rbo_unrank(rop, leftSum, left, d, leftRank, start, ctx);
  inner_rbo_unrank(rop, rightSum, right, d, rightRank, start + left, ctx);
}

void rbo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  inner_rbo_unrank(rop, n, k, d, r, 0, ctx);
}

uintx rbo_rank(const uint32_t n, const uint32_t k, const uint32_t d,
               const uint32_t *comb, const bic_ctx_t ctx) {
  if (k == 1) {
    return 0;
  }

  uint32_t left = (uint32_t)(k / 2);
  uint32_t right = k - left;

  const uint32_t *xl = comb;
  uint32_t leftSum = 0;
  for (uint32_t i = 0; i < left; ++i) {
    leftSum += xl[i];
  }

  const uint32_t *xr = comb + left;
  uint32_t rightSum = 0;
  for (uint32_t i = 0; i < right; ++i) {
    rightSum += xr[i];
  }

  uintx case3 = 0;
  for (uint32_t s = 0; s < leftSum; ++s) {
    case3 += ctx->comp(s, left, d, ctx) * ctx->comp(n - s, right, d, ctx);
  }

  uintx case5 =
      rbo_rank(leftSum, left, d, xl, ctx) * ctx->comp(rightSum, right, d, ctx);
  uintx case7 = rbo_rank(rightSum, right, d, xr, ctx);

  return case3 + case5 + case7;
}

static void topo_balanced(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k / 2;
  *k_r = k - *k_l;
}

static void iter_rbo_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                          uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->current_i = min_wl;
  s->val_1 = max_wl;
}

static int iter_rbo_next(bic_iter_state_t *s, uint32_t *w) {
  if (s->current_i > s->val_1)
    return 0;
  *w = s->current_i++;
  return 1;
}

static const bic_iter_strategy_t rbo_strat = {
    .name = "RBO",
    .topo = topo_balanced,
    .iter_init = iter_rbo_init,
    .iter_next = iter_rbo_next,
    .map = map_std,
    .unmap = unmap_std,
};

const bic_iter_strategy_t *bic_strat_rbo(void) { return &rbo_strat; }

void rbo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                        const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_rbo();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx rbo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                       const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_rbo();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}

void inner_rbo_mo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const uint32_t start,
                         const bic_ctx_t ctx) {
  uintx rank = r;

  if (k == 1) {
    rop[start] = n;
    return;
  }

  uint32_t left = (uint32_t)(k / 2);
  uint32_t right = k - left;

  uint32_t min_wl, max_wl;
  compute_bounds(n, left, right, d, &min_wl, &max_wl);

  uint32_t mid = (min_wl + max_wl) / 2;
  uint32_t leftSum = mid;

  uint32_t count_steps = 0;
  while (1) {
    uint32_t jump = (count_steps + 1) / 2;
    int sign = count_steps % 2;

    if (sign) {
      if (jump > mid) {
        leftSum = min_wl - 1;
      } else {
        leftSum = mid - jump;
      }
    } else {
      leftSum = mid + jump;
    }
    count_steps++;

    if (leftSum < min_wl || leftSum > max_wl) {
      if (jump > (max_wl > mid ? max_wl - mid : 0) &&
          jump > (mid > min_wl ? mid - min_wl : 0)) {
        break;
      }
      continue;
    }

    uint32_t rightSum = n - leftSum;
    uintx rightPoints = ctx->comp(rightSum, right, d, ctx);
    uintx count = ctx->comp(leftSum, left, d, ctx) * rightPoints;

    if (rank < count) {
      uintx leftRank = rank / rightPoints;
      uintx rightRank = rank % rightPoints;
      inner_rbo_mo_unrank(rop, leftSum, left, d, leftRank, start, ctx);
      inner_rbo_mo_unrank(rop, rightSum, right, d, rightRank, start + left,
                          ctx);
      return;
    }
    rank -= count;
  }
}

void rbo_mo_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                   const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  inner_rbo_mo_unrank(rop, n, k, d, r, 0, ctx);
}

uintx rbo_mo_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t *comb, const bic_ctx_t ctx) {
  if (k == 1) {
    return 0;
  }

  uint32_t left = (uint32_t)(k / 2);
  uint32_t right = k - left;

  const uint32_t *xl = comb;
  uint32_t leftSum = 0;
  for (uint32_t i = 0; i < left; ++i) {
    leftSum += xl[i];
  }

  const uint32_t *xr = comb + left;
  uint32_t rightSum = n - leftSum;

  uint32_t min_wl, max_wl;
  compute_bounds(n, left, right, d, &min_wl, &max_wl);
  uint32_t mid = (min_wl + max_wl) / 2;

  uintx ranks_before = 0;
  uint32_t count_steps = 0;

  while (1) {
    uint32_t val;
    uint32_t jump = (count_steps + 1) / 2;
    int sign = count_steps % 2;

    if (sign) {
      if (jump > mid)
        val = min_wl - 1;
      else
        val = mid - jump;
    } else {
      val = mid + jump;
    }
    count_steps++;

    if (val < min_wl || val > max_wl) {
      if (jump > (max_wl > mid ? max_wl - mid : 0) &&
          jump > (mid > min_wl ? mid - min_wl : 0)) {
        break;
      }
      continue;
    }

    if (val == leftSum) {
      uintx c_l = rbo_mo_rank(leftSum, left, d, xl, ctx);
      uintx c_r = rbo_mo_rank(rightSum, right, d, xr, ctx);
      uintx vol_r = ctx->comp(rightSum, right, d, ctx);
      return ranks_before + c_l * vol_r + c_r;
    }

    ranks_before +=
        ctx->comp(val, left, d, ctx) * ctx->comp(n - val, right, d, ctx);
  }
  return 0;
}

static void iter_rbo_mo_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                             uint32_t k_r, uint32_t d) {
  uint32_t min_wl, max_wl;
  compute_bounds(n, k_l, k_r, d, &min_wl, &max_wl);

  s->n = min_wl;
  s->val_1 = max_wl;
  s->val_2 = (min_wl + max_wl) / 2;
  s->current_i = 0;
}

static int iter_rbo_mo_next(bic_iter_state_t *s, uint32_t *w) {
  uint32_t min_wl = s->n;
  uint32_t max_wl = s->val_1;
  uint32_t mid = s->val_2;

  while (1) {
    uint32_t jump = (s->current_i + 1) / 2;
    int sign = s->current_i % 2;
    uint32_t val;

    if (sign) {
      if (jump > mid)
        val = min_wl - 1;
      else
        val = mid - jump;
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

static const bic_iter_strategy_t rbo_mo_strat = {
    .name = "RBO-MO",
    .topo = topo_balanced,
    .iter_init = iter_rbo_mo_init,
    .iter_next = iter_rbo_mo_next,
    .map = map_std,
    .unmap = unmap_std,
};

const bic_iter_strategy_t *bic_strat_rbo_mo(void) { return &rbo_mo_strat; }

void rbo_mo_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                           const uint32_t d, const uintx r,
                           const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_rbo_mo();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx rbo_mo_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                          const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_rbo_mo();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
