#include "gray.h"
#include "math.h"

#include "types.h"

void gray_unrank(uint32_t *rop, const uint16_t n, const uint16_t k,
                 const uint16_t d, const uintx r, const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = r;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0; count = ctx->comp(it_n - part, i, d, ctx), rank >= count;
         ++part, rank -= count) {
    }
    if (part & 1U) {
      rank = count - 1 - rank;
    }
  }

  rop[0] = it_n;
}

uintx gray_rank(const uint16_t n, const uint16_t k, const uint16_t d,
                const uint32_t *comb, const bic_ctx_t ctx) {
  uint16_t it_n = n;
  uintx rank = 0;
  uint16_t part = 0;

  for (uint16_t i = k - 1, p = n & 1U, parity = 0;
       part = comb[i], parity = it_n & 1U, i > 0; --i, it_n -= part) {
    if (parity == p) {
      for (uint16_t j = 0; j < part; ++j) {
        rank += ctx->comp(it_n - j, i, d, ctx);
      }
    } else {
      for (uint16_t j = part + 1; j < min(it_n, d) + 1; ++j) {
        rank += ctx->comp(it_n - j, i, d, ctx);
      }
    }
  }

  return rank;
}
