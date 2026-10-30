#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "api.h"
#include "cache/sparse_acc.h"

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_precompute(65536, 1, 1, ctx) == 0);
  assert(ctx->comb_cache->rows == 65537);
  assert(ctx->bin(65536, 1, ctx) == 65536);
  assert(ctx->comp(0, 0, 1, ctx) == 1);
  assert(ctx->comp(65536, 1, 1, ctx) == 0);
  bic_free_precomputed(ctx);

  cache_t *cache = (cache_t *)calloc(1, sizeof(*cache));
  assert(cache != NULL);
  cache->rows = 65537;
  cache->cols = 1;
  cache->meta = (bic_meta_t *)calloc(cache->rows, sizeof(*cache->meta));
  assert(cache->meta != NULL);
  cache->meta[65536].base = 70000;
  cache->meta[65536].length = 3;
  ctx->small_acc_cache = cache;
  uint32_t lower, upper;
  assert(small_acc_get_bounds(65536, 0, &lower, &upper, ctx) == 0);
  assert(lower == 70000 && upper == 70003);
  bic_ctx_destroy(ctx);

  puts("cache bounds regressions passed");
  return 0;
}
