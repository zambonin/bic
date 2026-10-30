#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "generic.h"
#include "generic/rbo.h"
#include "order/rbo.h"

static bic_math_comp_func_t original_comp;
static unsigned long calls;

static uintx feasible_count(uint32_t n, uint32_t k, uint32_t d,
                            const bic_ctx_t ctx) {
  assert((uint64_t)n <= (uint64_t)k * d);
  ++calls;
  return original_comp(n, k, d, ctx);
}

int main(void) {
  unsigned long checked = 0;
  unsigned long specialized_calls = 0, generic_calls = 0;
  for (uint32_t d = 0; d <= 3; ++d) {
    uint32_t counts[8][22] = {{0}};
    counts[0][0] = 1;
    for (uint32_t k = 1; k <= 7; ++k) {
      for (uint32_t n = 0; n <= k * d; ++n) {
        for (uint32_t a = 0; a <= d && a <= n; ++a) {
          counts[k][n] += counts[k - 1][n - a];
        }
      }
      for (uint32_t n = 0; n <= k * d; ++n) {
        const uint32_t total = counts[k][n];
        if (total > 1200 || !bic_params_supported(n, k, d, BIC_CACHE_NONE)) {
          continue;
        }
        bic_ctx_t ctx = bic_ctx_init();
        assert(ctx != NULL);
        assert(bic_ctx_set_order(BIC_ORDER_RBO, ctx) == 0);
        assert(bic_precompute(n, k, d, ctx) == 0);
        assert(ctx->comp(n, k, d, ctx) == total);
        original_comp = ctx->comp;
        ctx->comp = feasible_count;
        bool seen[16384] = {false};
        for (uint32_t r = 0; r < total; ++r) {
          uint32_t a[7], b[7];
          calls = 0;
          rbo_unrank(a, n, k, d, r, ctx);
          specialized_calls += calls;
          calls = 0;
          rbo_generic_unrank(b, n, k, d, r, ctx);
          generic_calls += calls;
          assert(memcmp(a, b, k * sizeof(*a)) == 0);
          uint32_t sum = 0, key = 0;
          for (uint32_t i = 0; i < k; ++i) {
            assert(a[i] <= d);
            sum += a[i];
            key = key * (d + 1) + a[i];
          }
          assert(sum == n && !seen[key]);
          seen[key] = true;
          assert(rbo_rank(n, k, d, a, ctx) == r);
          assert(rbo_generic_rank(n, k, d, b, ctx) == r);
          ++checked;
        }
        bic_ctx_destroy(ctx);
      }
    }
  }
  uint32_t lo, hi;
  compute_bounds(UINT32_MAX, UINT32_MAX / 2, UINT32_MAX / 2 + 1, UINT32_MAX,
                 &lo, &hi);
  assert(lo == 0 && hi == UINT32_MAX);
  compute_bounds(20, 5, 5, 2, &lo, &hi);
  assert(lo == 10 && hi == 10);
  printf("RBO: %lu ranks; unrank count calls specialized=%lu generic=%lu\n",
         checked, specialized_calls, generic_calls);
  return 0;
}
