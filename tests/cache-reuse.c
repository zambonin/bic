#include <assert.h>
#include <stdio.h>

#include "api.h"

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);

  assert(bic_precompute(2, 2, 2, ctx) == 0);
  assert(ctx->comp(2, 2, 2, ctx) == 3);
  assert(bic_precompute(5, 2, 5, ctx) == 0);
  assert(ctx->bin(5, 2, ctx) == 10);
  assert(ctx->comp(5, 2, 5, ctx) == 6);
  assert(bic_precompute(5, 2, 1, ctx) == 0);
  assert(ctx->comp(5, 2, 1, ctx) == 0);

  assert(bic_ctx_set_cache(BIC_CACHE_NONE, ctx) == 0);
  assert(ctx->bin_cache == NULL && ctx->comb_cache == NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_precompute(2, 2, 2, ctx) == 0);
  cache_t *dense = ctx->comb_cache;
  assert(bic_ctx_set_order(BIC_ORDER_GRAY, ctx) == 0);
  assert(ctx->comb_cache == dense && ctx->prepared);
  for (uint32_t rank = 0; rank < 3; ++rank) {
    uint32_t parts[2];
    bic_unrank(parts, 2, 2, 2, rank, ctx);
    assert(parts[0] + parts[1] == 2);
    assert(bic_rank(2, 2, 2, parts, ctx) == rank);
  }

  bic_free_precomputed(ctx);
  assert(bic_ctx_set_order(BIC_ORDER_RBO, ctx) == 0);
  assert(bic_precompute(6, 6, 2, ctx) == 0);
  assert(ctx->comb_cache->col_map != NULL);
  assert(bic_ctx_set_order(BIC_ORDER_GRAY, ctx) == 0);
  assert(ctx->bin_cache == NULL && ctx->comb_cache == NULL);
  assert(!ctx->prepared);

  assert(bic_ctx_set_cache(BIC_CACHE_SMALL_ACC, ctx) == 0);
  assert(bic_precompute(3, 2, 2, ctx) == 0);
  bic_free_precomputed(ctx);
  assert(ctx->dir == compute_dir);

  bic_ctx_destroy(ctx);
  puts("cache reuse regression passed");
  return 0;
}
