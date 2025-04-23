#include "rbo.h"

void inner_rbo(uint32_t *rop, const uint16_t n, const uint16_t k,
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

  inner_rbo(rop, leftSum, left, d, leftRank, start);
  inner_rbo(rop, rightSum, right, d, rightRank, start + left);
}

void rbo(uint32_t *rop, const uint16_t n, const uint16_t k, const uint16_t d,
         const uintx r) {
  inner_rbo(rop, n, k, d, r, 0);
}
