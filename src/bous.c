#include "bous.h"
#include "generic.h"
#include "math.h"
#include "types.h"

void inner_bous_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
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

  for (uint32_t j = 0; j <= n / 2; ++j) {
    uint32_t S[2] = {j, (uint32_t)(n - j)};
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

void bous_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                 const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  inner_bous_unrank(rop, n, k, d, r, 0, ctx);
}

uintx bous_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                const uint32_t *comb, const bic_ctx_t ctx) {
  if (k == 1) {
    return 0;
  }

  uint32_t left = (uint32_t)(k / 2);
  uint32_t right = k - left;

  uint32_t leftSum = 0;
  for (uint32_t i = 0; i < left; ++i) {
    leftSum += comb[i];
  }

  uint32_t rightSum = n - leftSum;

  uintx c_pre = 0;
  bool found = false;

  for (uint32_t j = 0; j <= n / 2; ++j) {
    uint32_t S[2] = {j, (uint32_t)(n - j)};
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

static void topo_balanced(uint32_t k, uint32_t *k_l, uint32_t *k_r) {
  *k_l = k / 2;
  *k_r = k - *k_l;
}

static void iter_bous_init(bic_iter_state_t *s, uint32_t n, uint32_t k_l,
                           uint32_t k_r, uint32_t d) {
  (void)k_l;
  (void)k_r;
  (void)d;

  s->n = n;
  s->val_1 = 0;
  s->val_2 = n / 2;
  s->current_i = 0;
}

static int iter_bous_next(bic_iter_state_t *s, uint32_t *w) {
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

static const bic_iter_strategy_t bous_strat = {
    .name = "Bous",
    .topo = topo_balanced,
    .iter_init = iter_bous_init,
    .iter_next = iter_bous_next,
    .map = map_std,
    .unmap = unmap_std,
};

const bic_iter_strategy_t *bic_strat_bous(void) { return &bous_strat; }

void bous_generic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                         const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_bous();
  bic_generic_unrank(rop, n, k, d, &r, strat, ctx);
}

uintx bous_generic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *comb, const bic_ctx_t ctx) {
  const bic_iter_strategy_t *strat = bic_strat_bous();
  return bic_generic_rank(n, k, d, comb, strat, ctx);
}
