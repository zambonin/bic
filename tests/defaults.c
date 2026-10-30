#include <assert.h>
#include <stdio.h>

#include "api.h"

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  bic_ctx_t fresh = bic_ctx_init();
  assert(ctx != NULL);
  assert(fresh != NULL);
  assert(bic_ctx_set_order(BIC_ORDER_GRAY, ctx) == 0);
  assert(bic_ctx_set_defaults(ctx) == 0);
  assert(ctx->order == BIC_ORDER_COLEX);
  assert(ctx->unrank_alg == BIC_ALG_DEFAULT);
  assert(ctx->rank_alg == BIC_RANK_ALG_DEFAULT);
  assert(ctx->unrank == fresh->unrank);
  assert(ctx->rank == fresh->rank);
  bic_ctx_destroy(ctx);
  bic_ctx_destroy(fresh);
  puts("context defaults regression passed");
  return 0;
}
