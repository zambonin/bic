#include "order/rbo.h"
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

  const uint64_t capacity = (uint64_t)left * d;
  const uint32_t limit = capacity < n ? (uint32_t)capacity : n;
  for (uintx count = 0; leftSum <= limit; ++leftSum, rank -= count) {
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
