#include "rbo.h"
#include "math.h"

#include "types.h"

void inner_rbo_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                      const uint16_t k, const uint16_t d, const uintx r,
                      const uint16_t start) {
  uintx rank = r;

  if (k == 1) {
    rop[start] = n;
    return;
  }

  uint16_t left = (uint16_t)(k / 2);
  uint16_t right = k - left;

  uint16_t leftSum = 0;
  uint16_t rightSum = 0;
  uintx rightPoints = 0;

  for (uintx count = 0; leftSum <= min(n, left * d); ++leftSum, rank -= count) {
    rightSum = n - leftSum;
    rightPoints = ctx->comp(ctx, rightSum, right, d);
    count = ctx->comp(ctx, leftSum, left, d) * rightPoints;
    if (rank < count) {
      break;
    }
  }

  uintx leftRank = rank / rightPoints;
  uintx rightRank = rank % rightPoints;

  inner_rbo_unrank(ctx, rop, leftSum, left, d, leftRank, start);
  inner_rbo_unrank(ctx, rop, rightSum, right, d, rightRank, start + left);
}

void rbo_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                const uint16_t k, const uint16_t d, const uintx r) {
  inner_rbo_unrank(ctx, rop, n, k, d, r, 0);
}

uintx rbo_rank(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
               const uint16_t d, const uint32_t *comb) {
  if (k == 1) {
    return 0;
  }

  uint16_t left = (uint16_t)(k / 2);
  uint16_t right = k - left;

  const uint32_t *xl = comb;
  uint16_t leftSum = 0;
  for (int i = 0; i < left; ++i) {
    leftSum += xl[i];
  }

  const uint32_t *xr = comb + left;
  uint16_t rightSum = 0;
  for (int i = 0; i < right; ++i) {
    rightSum += xr[i];
  }

  uintx case3 = 0;
  for (uint16_t s = 0; s < leftSum; ++s) {
    case3 += ctx->comp(ctx, s, left, d) * ctx->comp(ctx, n - s, right, d);
  }

  uintx case5 =
      rbo_rank(ctx, leftSum, left, d, xl) * ctx->comp(ctx, rightSum, right, d);
  uintx case7 = rbo_rank(ctx, rightSum, right, d, xr);

  return case3 + case5 + case7;
}
