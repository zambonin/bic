#include <assert.h>
#include <stdio.h>

#include "api.h"
#include "cache/comb.h"
#include "cache/common.h"

int main(void) {
  if (SIZE_MAX > UINT32_MAX) {
    assert(comb_cache_length(65536, 65536) == (size_t)4294967296ULL);
  }
  assert(!cache_product_fits(SIZE_MAX, 2, sizeof(uintx)));

  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_precompute(UINT32_MAX - 1, 1, 1, ctx) != 0);
  assert(ctx->bin_cache == NULL && ctx->comb_cache == NULL);
  assert(bic_precompute(4, 2, 2, ctx) == 0);
  bic_ctx_destroy(ctx);

  puts("cache allocation regressions passed");
  return 0;
}
