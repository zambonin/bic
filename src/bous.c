#include "bous.h"
#include "math.h"
#include "types.h"

void inner_bous_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                       const uint16_t d, const uintx r, const uint16_t start,
                       const bic_ctx_t ctx) {
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

  for (uint16_t j = 0; j <= n / 2; ++j) {
    uint16_t S[2] = {j, (uint16_t)(n - j)};
    uint8_t S_len = (j == n - j) ? 1 : 2;

    for (uint8_t i = 0; i < S_len; ++i) {
      leftSum = S[i];
      rightSum = n - leftSum;
      rightPoints = ctx->comp(rightSum, right, d, ctx);
      uintx count = ctx->comp(leftSum, left, d, ctx) * rightPoints;

      if (rank < count) {
        inner_bous_unrank(rop, leftSum, left, d, rank / rightPoints, start,
                          ctx);
        inner_bous_unrank(rop, rightSum, right, d, rank % rightPoints,
                          start + left, ctx);
        return;
      }
      rank -= count;
    }
  }
}

void bous_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                 const uint16_t d, const uintx r, const bic_ctx_t ctx) {
  inner_bous_unrank(rop, n, k, d, r, 0, ctx);
}

uintx bous_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                const uint32_t *comb, const bic_ctx_t ctx) {
  if (k == 1) {
    return 0;
  }

  uint16_t left = (uint16_t)(k / 2);
  uint16_t right = k - left;

  uint16_t leftSum = 0;
  for (int i = 0; i < left; ++i) {
    leftSum += comb[i];
  }

  uint16_t rightSum = n - leftSum;

  uintx c_pre = 0;
  bool found = false;

  for (uint16_t j = 0; j <= n / 2; ++j) {
    uint16_t S[2] = {j, (uint16_t)(n - j)};
    uint8_t S_len = (j == n - j) ? 1 : 2;

    for (uint8_t i = 0; i < S_len; ++i) {
      if (S[i] == leftSum) {
        found = true;
        break;
      }
      c_pre +=
          ctx->comp(S[i], left, d, ctx) * ctx->comp(n - S[i], right, d, ctx);
    }
    if (found) {
      break;
    }
  }

  uintx c_l = bous_rank(leftSum, left, d, comb, ctx) *
              ctx->comp(rightSum, right, d, ctx);
  uintx c_r = bous_rank(rightSum, right, d, comb + left, ctx);

  return c_pre + c_l + c_r;
}
