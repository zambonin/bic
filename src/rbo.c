#include "rbo.h"
#include "math.h"
#include "utils.h"

void inner_rbo_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                      const uint16_t d, const uintx r, const uint16_t start) {
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
    rightPoints = bic(rightSum, right, d);
    count = bic(leftSum, left, d) * rightPoints;
    if (rank < count) {
      break;
    }
  }

  uintx leftRank = rank / rightPoints;
  uintx rightRank = rank % rightPoints;

  inner_rbo_unrank(rop, leftSum, left, d, leftRank, start);
  inner_rbo_unrank(rop, rightSum, right, d, rightRank, start + left);
}

void rbo_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                const uint16_t d, const uintx r) {
  inner_rbo_unrank(rop, n, k, d, r, 0);
}

uintx rbo_rank(const uint16_t n, const uint16_t k, const uint16_t d,
               const uint32_t *comb) {
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
    case3 += bic(s, left, d) * bic(n - s, right, d);
  }

  uintx case5 = rbo_rank(leftSum, left, d, xl) * bic(rightSum, right, d);
  uintx case7 = rbo_rank(rightSum, right, d, xr);

  return case3 + case5 + case7;
}
