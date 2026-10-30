#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "api.h"
#include "math.h"

static uint64_t count_compositions(uint32_t n, uint32_t k, uint32_t d) {
  uint64_t dp[13] = {1};
  uint64_t next[13];
  for (uint32_t part = 0; part < k; ++part) {
    for (uint32_t sum = 0; sum <= k * d; ++sum) {
      next[sum] = 0;
      for (uint32_t value = 0; value <= d && value <= sum; ++value) {
        next[sum] += dp[sum - value];
      }
    }
    for (uint32_t sum = 0; sum <= k * d; ++sum) {
      dp[sum] = next[sum];
    }
  }
  return dp[n];
}

int main(void) {
  for (uint32_t k = 1; k <= 4; ++k) {
    for (uint32_t d = 1; d <= 3; ++d) {
      for (uint32_t n = 0; n <= k * d; ++n) {
        const uint64_t count = count_compositions(n, k, d);
        assert(compute_bic_no_cache(n, k, d) == count);
        for (uint8_t order = 0; order < BIC_ORDER_LENGTH; ++order) {
          for (uint8_t cache = BIC_CACHE_NONE; cache < BIC_CACHE_LENGTH;
               ++cache) {
            bic_ctx_t ctx = bic_ctx_init();
            assert(ctx != NULL);
            assert(bic_ctx_set_order(order, ctx) == 0);
            assert(bic_ctx_set_cache(cache, ctx) == 0);
            assert(bic_precompute(n, k, d, ctx) == 0);
            bool seen[256] = {false};
            for (uint64_t rank = 0; rank < count; ++rank) {
              uint32_t comp[4] = {0};
              bic_unrank(comp, n, k, d, rank, ctx);
              uint32_t sum = 0;
              uint32_t code = 0;
              uint32_t place = 1;
              for (uint32_t i = 0; i < k; ++i) {
                assert(comp[i] <= d);
                sum += comp[i];
                code += comp[i] * place;
                place *= d + 1;
              }
              assert(sum == n && code < 256 && !seen[code]);
              seen[code] = true;
              assert(bic_rank(n, k, d, comp, ctx) == rank);
            }
            bic_ctx_destroy(ctx);
          }
        }
      }
    }
  }
  const uint32_t wide_d = UINT32_C(1073741824);
  assert(bic_params_supported(1, 8, wide_d, BIC_CACHE_NONE));
  assert(compute_bic_no_cache(1, 8, wide_d) == 8);
  for (uint8_t order = 0; order < BIC_ORDER_LENGTH; ++order) {
    bic_ctx_t ctx = bic_ctx_init();
    assert(ctx != NULL);
    assert(bic_ctx_set_order(order, ctx) == 0);
    assert(bic_precompute(1, 8, wide_d, ctx) == 0);
    unsigned seen = 0;
    for (uint32_t rank = 0; rank < 8; ++rank) {
      uint32_t comp[8] = {0};
      bic_unrank(comp, 1, 8, wide_d, rank, ctx);
      uint32_t sum = 0;
      unsigned code = 0;
      for (uint32_t i = 0; i < 8; ++i) {
        assert(comp[i] <= 1);
        sum += comp[i];
        code |= comp[i] << i;
      }
      assert(sum == 1 && !(seen & code));
      seen |= code;
      assert(bic_rank(1, 8, wide_d, comp, ctx) == rank);
    }
    assert(seen == 255);
    bic_ctx_destroy(ctx);
  }
  puts("independent composition oracle passed");
  return 0;
}
