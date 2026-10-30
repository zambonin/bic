#include <assert.h>
#include <errno.h>
#include <stdio.h>

#include "api.h"
#include "strat.h"

static void assert_impossible_strategies_clear_outputs(const uint32_t m,
                                                       const uint32_t k) {
  uint32_t n = 123;
  uint32_t d = 456;
  mingen(m, &n, k, &d);
  assert(n == UINT32_MAX && d == UINT32_MAX && errno == ERANGE);

  n = 123;
  d = 456;
  minver(m, &n, k, &d);
  assert(n == UINT32_MAX && d == UINT32_MAX && errno == ERANGE);
}

int main(void) {
  assert_impossible_strategies_clear_outputs(1, 1);
  assert_impossible_strategies_clear_outputs(32, 2);

  assert(!bic_params_supported(1, 2, UINT32_MAX, BIC_CACHE_NONE));
  assert(!bic_params_supported(UINT32_MAX - 1, 2, 1, BIC_CACHE_NONE));
  assert(bic_params_supported(1, 10007, 1, BIC_CACHE_NONE));
  assert(compute_bic_no_cache(1, 10007, 1) == 10007);

  errno = 0;
  assert(compute_bic_no_cache(1, 2, UINT32_MAX) == 0);
  assert(errno == ERANGE);

#if defined(BITINT) || defined(BOOST_FIX_INT)
  assert(!bic_params_supported(1000, 500, 1, BIC_CACHE_NONE));
  assert(!bic_params_supported(1, 10007, 1, BIC_CACHE_BIN));
  errno = 0;
  assert(compute_bin(1000, 500, NULL) == 0);
  assert(errno == ERANGE);
#endif

  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_precompute(1, 2, UINT32_MAX, ctx) != 0);
  assert(ctx->bin_cache == NULL && ctx->comb_cache == NULL);
  assert(bic_precompute(3, 2, 2, ctx) == 0);
  bic_ctx_destroy(ctx);

  puts("parameter domain regressions passed");
  return 0;
}
