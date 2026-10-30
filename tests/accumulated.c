#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "cache/common.h"
#include "search.h"

int main(void) {
  const uintx repeated[] = {0, 0, 2, 2, 5};
  uintx key = 0;
  assert(uintx_upper_bound(NULL, 0, &key) == 0);
  assert(uintx_upper_bound(repeated, 5, &key) == 2);
  key = 2;
  assert(uintx_upper_bound(repeated, 5, &key) == 4);
  key = 4;
  assert(uintx_upper_bound(repeated, 5, &key) == 4);
  key = 5;
  assert(uintx_upper_bound(repeated, 5, &key) == 5);

  unsigned long checked = 0;
  for (uint32_t d = 0; d <= 3; ++d) {
    uint32_t counts[7][19] = {{0}};
    counts[0][0] = 1;
    uint32_t tuples = 1;
    for (uint32_t k = 1; k <= 6; ++k) {
      tuples *= d + 1;
      for (uint32_t n = 0; n <= k * d; ++n) {
        for (uint32_t a = 0; a <= d && a <= n; ++a) {
          counts[k][n] += counts[k - 1][n - a];
        }
      }
      for (uint32_t n = 0; n <= k * d; ++n) {
        for (uint8_t cache = 0; cache < BIC_CACHE_LENGTH; ++cache) {
          bic_ctx_t ctx = bic_ctx_init();
          assert(ctx != NULL);
          assert(bic_ctx_set_cache(cache, ctx) == 0);
          if (bic_precompute(n, k, d, ctx)) {
            bic_ctx_destroy(ctx);
            continue;
          }
          if (ctx->acc_cache != NULL) {
            const cache_t *c = ctx->acc_cache;
            size_t entries = 0;
            for (uint32_t row = 0; row <= n; ++row) {
              for (uint32_t col = 0; col < k; ++col) {
                uint32_t len;
                const uintx *slice = acc_cache_get_ptr(row, col, &len, ctx);
                assert(slice != NULL && len == min(row, d) + 1);
                uint32_t sum = 0;
                assert(ctx->dir(row, col, d, 0, ctx) == 0);
                for (uint32_t a = 0; a < len; ++a) {
                  sum += counts[col][row - a];
                  assert(slice[a] == sum);
                  assert(ctx->dir(row, col, d, a + 1, ctx) == sum);
                }
                entries += len;
              }
            }
            assert(c->length == entries);
            uint32_t len = 123;
            assert(acc_cache_get_ptr(n + 1, k - 1, &len, ctx) == NULL);
            assert(len == 0);
          }
          for (uint8_t alg = 0; alg < BIC_UNRANK_ALG_LENGTH; ++alg) {
            if (bic_ctx_set_unrank_alg(alg, ctx)) {
              continue;
            }
            uint32_t r = 0;

            for (uint32_t word = 0; word < tuples; ++word) {
              uint32_t expected[6], actual[6], q = word, sum = 0;
              for (uint32_t i = 0; i < k; ++i) {
                expected[i] = q % (d + 1);
                q /= d + 1;
                sum += expected[i];
              }
              if (sum != n) {
                continue;
              }
              bic_unrank(actual, n, k, d, r, ctx);
              assert(memcmp(expected, actual, k * sizeof(*actual)) == 0);
              assert(bic_rank(n, k, d, actual, ctx) == r);
              if (bic_ctx_set_rank_alg(BIC_RANK_ALG_AD, ctx) == 0) {
                assert(bic_rank(n, k, d, actual, ctx) == r);
                assert(bic_ctx_set_rank_alg(BIC_RANK_ALG_DEFAULT, ctx) == 0);
              }
              ++r;
              ++checked;
            }
            assert(r == counts[k][n]);
          }
          bic_ctx_destroy(ctx);
        }
      }
    }
  }

  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_order(BIC_ORDER_RBO, ctx) == 0);
  assert(bic_ctx_set_cache(BIC_CACHE_ACC, ctx) == 0);
  assert(bic_precompute(0, 1, 0, ctx) == 0);
  uint32_t len = 123;
  assert(acc_cache_get_ptr(0, 0, &len, ctx) == NULL && len == 0);
  bic_ctx_destroy(ctx);
  printf("Accumulated prefixes: %lu colex ranks, all algorithms/caches\n",
         checked);
  return 0;
}
