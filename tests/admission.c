#include "api.h"
#include "cache/common.h"
#include <assert.h>
#include <errno.h>
#include <stdio.h>

static unsigned exact_count(unsigned n, unsigned k, unsigned d) {
  if (k == 0) {
    return n == 0;
  }
  unsigned count = 0;
  for (unsigned x = 0; x <= d && x <= n; ++x) {
    count += exact_count(n - x, k - 1, d);
  }
  return count;
}

static void check_class(unsigned n, unsigned k, unsigned d, bic_ctx_t ctx) {
  const unsigned count = exact_count(n, k, d);
  assert(ctx->comp(n, k, d, ctx) == count);
  bool seen[243] = {false};
  for (unsigned r = 0; r < count; ++r) {
    uint32_t parts[5] = {0};
    bic_unrank(parts, n, k, d, r, ctx);
    unsigned sum = 0, code = 0;
    for (unsigned i = 0; i < k; ++i) {
      assert(parts[i] <= d);
      sum += parts[i];
      code = code * (d + 1) + parts[i];
    }
    assert(sum == n && !seen[code]);
    seen[code] = true;
    assert(bic_rank(n, k, d, parts, ctx) == r);
  }
}

int main(void) {
  assert(cache_band_endpoint(-1e300, INT32_MAX, false) == 0);
  assert(cache_band_endpoint(1e300, INT32_MAX, true) == INT32_MAX);
  assert(cache_band_endpoint(3.75, 10, false) == 3);
  assert(cache_band_endpoint(__builtin_nan(""), 10, false) == 0);
  assert(cache_band_endpoint(__builtin_nan(""), 10, true) == 10);
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL && !ctx->prepared);
  uint32_t parts[2] = {99, 99};
  errno = 0;
  bic_unrank(parts, 1, 2, 1, 0, ctx);
  assert(errno == EINVAL && parts[0] == 99 && parts[1] == 99);
  assert(bic_precompute(1, 2, 1, ctx) == 0);
  assert(ctx->prepared && !ctx->precomputed);
  assert(bic_ctx_set_unrank_alg(BIC_ALG_PS, ctx) == 0);
  check_class(1, 2, 1, ctx);
  assert(bic_precompute(UINT32_MAX, 2, 1, ctx) != 0);
  assert(ctx->prepared && ctx->precomputed_n == 1);
  check_class(1, 2, 1, ctx);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(!ctx->prepared);
  assert(bic_precompute(1, 2, 1, ctx) == 0);
  assert(ctx->prepared && ctx->precomputed);
  assert(bic_ctx_set_order(BIC_ORDER_RBO, ctx) == 0);
  assert(ctx->prepared);
  check_class(1, 2, 1, ctx);
  bic_free_precomputed(ctx);
  assert(!ctx->prepared && !ctx->precomputed);
  assert(bic_precompute(0, 0, 0, ctx) == 0);
  assert(ctx->prepared && !ctx->precomputed);
  assert(bic_rank(0, 0, 0, NULL, ctx) == 0);
  bic_ctx_destroy(ctx);

#if defined(BITINT) && BITINT == 16
  ctx = bic_ctx_init();
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_precompute(20, 10, 2, ctx) == 0);
  for (uint8_t alg = BIC_ALG_PS; alg <= BIC_ALG_AD; ++alg) {
    assert(bic_ctx_set_unrank_alg(alg, ctx) != 0);
    assert(ctx->unrank_alg == BIC_ALG_DEFAULT && ctx->prepared);
  }
  assert(bic_ctx_set_rank_alg(BIC_RANK_ALG_AD, ctx) != 0);
  assert(ctx->rank_alg == BIC_RANK_ALG_DEFAULT && ctx->prepared);
  uint32_t endpoint[10] = {0};
  bic_unrank(endpoint, 20, 10, 2, 0, ctx);
  for (unsigned i = 0; i < 10; ++i) {
    assert(endpoint[i] == 2);
  }
  assert(bic_ctx_set_unrank_alg(BIC_ALG_GENERIC, ctx) == 0);
  bic_unrank(endpoint, 20, 10, 2, 0, ctx);
  for (unsigned i = 0; i < 10; ++i) {
    assert(endpoint[i] == 2);
  }
  bic_ctx_destroy(ctx);
#endif

  unsigned classes = 0, rejected = 0;
  for (unsigned k = 1; k <= 5; ++k) {
    for (unsigned d = 0; d <= 2; ++d) {
      for (unsigned n = 0; n <= k * d; ++n) {
        for (uint8_t order = 0; order < BIC_ORDER_LENGTH; ++order) {
          for (uint8_t u = 0; u < BIC_UNRANK_ALG_LENGTH; ++u) {
            for (uint8_t r = 0; r < BIC_RANK_ALG_LENGTH; ++r) {
              for (uint8_t cache = 0; cache < BIC_CACHE_LENGTH; ++cache) {
                ctx = bic_ctx_init();
                assert(ctx != NULL);
                assert(bic_ctx_set_order(order, ctx) == 0);
                if (bic_ctx_set_unrank_alg(u, ctx) ||
                    bic_ctx_set_rank_alg(r, ctx)) {
                  bic_ctx_destroy(ctx);
                  continue;
                }
                assert(bic_ctx_set_cache(cache, ctx) == 0);
                if (bic_precompute(n, k, d, ctx)) {
                  ++rejected;
                } else {
                  check_class(n, k, d, ctx);
                  ++classes;
                }
                bic_ctx_destroy(ctx);
              }
            }
          }
        }
      }
    }
  }
  printf("admission: %u exact classes checked; %u configurations rejected\n",
         classes, rejected);
  return 0;
}
