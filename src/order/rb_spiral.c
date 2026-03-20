#include "order/rb_spiral.h"
#include "math.h"
#include "types.h"

static inline void compute_bounds(uint32_t n, uint32_t k_l, uint32_t k_r,
                                  uint32_t d, uint32_t *min_wl,
                                  uint32_t *max_wl) {
  uint64_t max_wr = (uint64_t)k_r * d;
  *min_wl = ((uint64_t)n > max_wr) ? n - (uint32_t)max_wr : 0;
  uint64_t max_wl_cap = (uint64_t)k_l * d;
  *max_wl = ((uint64_t)n < max_wl_cap) ? n : (uint32_t)max_wl_cap;
}

void inner_rb_spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                            const uint32_t d, const uintx r,
                            const uint32_t start, const bic_ctx_t ctx) {
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
      inner_rb_spiral_unrank(rop, leftSum, left, d, leftRank, start, ctx);
      inner_rb_spiral_unrank(rop, rightSum, right, d, rightRank, start + left,
                             ctx);
      return;
    }
    rank -= count;
  }
}

void rb_spiral_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                      const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  inner_rb_spiral_unrank(rop, n, k, d, r, 0, ctx);
}

uintx rb_spiral_rank(const uint32_t n, const uint32_t k, const uint32_t d,
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
      if (jump > mid) {
        val = min_wl - 1;
      } else {
        val = mid - jump;
      }
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
      uintx c_l = rb_spiral_rank(leftSum, left, d, xl, ctx);
      uintx c_r = rb_spiral_rank(rightSum, right, d, xr, ctx);
      uintx vol_r = ctx->comp(rightSum, right, d, ctx);
      return ranks_before + c_l * vol_r + c_r;
    }

    ranks_before +=
        ctx->comp(val, left, d, ctx) * ctx->comp(n - val, right, d, ctx);
  }

  return 0;
}
