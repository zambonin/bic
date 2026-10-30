#include <assert.h>
#include <stdio.h>

#include "api.h"

int main(void) {
  bic_ctx_destroy(NULL);
  for (uint8_t cache = BIC_CACHE_BIN; cache < BIC_CACHE_LENGTH; ++cache) {
    bic_ctx_t ctx = bic_ctx_init();
    assert(ctx != NULL);
    assert(bic_ctx_set_cache(cache, ctx) == 0);
    assert(bic_precompute(5, 3, 3, ctx) == 0);
    bic_ctx_destroy(ctx);
  }
  puts("context destruction regression passed");
  return 0;
}
