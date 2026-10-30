#include <assert.h>
#include <stdio.h>

#include "api.h"

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_BIN, ctx) == 0);
  assert(bic_precompute(2, 2, 2, ctx) == 0);
  assert(ctx->bin(2, 1, ctx) == 2);
  assert(ctx->bin(10007, 1, ctx) == 10007);
  assert(ctx->bin(5, 3, ctx) == 10);
  assert(ctx->bin(1, 2, ctx) == 0);
  bic_ctx_destroy(ctx);
  puts("binomial cache regression passed");
  return 0;
}
