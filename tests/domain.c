#include <assert.h>
#include <errno.h>
#include <stdint.h>

#include "api.h"
#include "arithmetic.h"
#include "extra.h"
#include "strat.h"

int main(void) {
  assert(bic_muldiv_exact(11480, 3, 42) == 820);
  assert(bic_signed_muldiv_exact(11480, 3, 42) == 820);
  assert(bic_signed_muldiv_exact(-11480, 3, 42) == -820);
  assert(bic_muldiv_exact(123, 0, 7) == 0);
  assert(bits_fit_bic(1, 0, 2) == 0 && errno == 0);
  assert(bits_fit_bic(2, 2, 2) == 2 && errno == 0);
  assert(bits_fit_bic(UINT32_MAX, 2, 2) == 0 && errno == ERANGE);
  assert(random_rank(1, 0, 2) == 0 && errno == EDOM);
  assert(random_rank(UINT32_MAX, 2, 2) == 0 && errno == ERANGE);
  assert(random_rank(0, 1, 2) == 0 && errno == 0);
  uintx value = 99;
  assert(bic_count_checked(&value, 1, 0, 2) && value == 0);
  assert(bic_compare_2_pow_m(0, 1, 0, 2) == BIC_COUNT_BELOW);
  assert(bic_compare_2_pow_m(0, 0, 0, 2) == BIC_COUNT_AT_LEAST);
  assert(!bic_count_checked(&value, UINT32_MAX, 2, 2) && value == 0);
  assert(errno == ERANGE);
  assert(bic_compare_2_pow_m(0, UINT32_MAX, 2, 2) == BIC_COUNT_UNSUPPORTED);
  assert(bic_count_checked(&value, 0, 2, 2) && value == 1 && errno == 0);
  uint32_t failed_n, failed_d;
  mingen(1, &failed_n, 1, &failed_d);
  assert(failed_n == UINT32_MAX && failed_d == UINT32_MAX && errno == ERANGE);
  minver(0, &failed_n, 1, &failed_d);
  assert(failed_n == 0 && failed_d == 0 && errno == 0);
#if !defined(BITINT) || BITINT >= 512

  const struct {
    uint32_t m, k, gen_n, gen_d, ver_n, ver_d;
  } wide_optima[] = {
      {260, 10, INT32_MAX, 686334074, 2061453506, 1877716767},
      {448, 17, INT32_MAX, 298129809, 1825486398, 1366537962},
  };
  for (size_t i = 0; i < sizeof(wide_optima) / sizeof(wide_optima[0]); ++i) {
    uint32_t n, d;
    assert(bic_params_supported(INT32_MAX, wide_optima[i].k, INT32_MAX,
                                BIC_CACHE_NONE));
    assert(bic_compare_2_pow_m(wide_optima[i].m, INT32_MAX, wide_optima[i].k,
                               INT32_MAX) == BIC_COUNT_AT_LEAST);
    mingen(wide_optima[i].m, &n, wide_optima[i].k, &d);
    assert(errno == 0 && n == wide_optima[i].gen_n &&
           d == wide_optima[i].gen_d);
    minver(wide_optima[i].m, &n, wide_optima[i].k, &d);
    assert(errno == 0 && n == wide_optima[i].ver_n &&
           d == wide_optima[i].ver_d);
  }
#endif
#if defined(BITINT) && BITINT == 16
  assert(bic_compare_2_pow_m(0, 10, 10, 2) == BIC_COUNT_UNSUPPORTED);
  assert(bic_compare_2_pow_m(0, 20, 10, 2) == BIC_COUNT_AT_LEAST);
  assert(bic_compare_2_pow_m(20, 3, 41, 3) == BIC_COUNT_BELOW);
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_params_supported(3, 41, 3, BIC_CACHE_NONE));
  assert(compute_bic_no_cache(3, 41, 3) == 12341);
  assert(bic_ctx_set_unrank_alg(BIC_ALG_PS, ctx) == 0);

  assert(bic_precompute(3, 41, 3, ctx) == 0);
  uint32_t parts[41] = {0};
  bic_unrank(parts, 3, 41, 3, 11480, ctx);
  uint32_t sum = 0;
  for (uint32_t i = 0; i < 41; ++i) {
    assert(parts[i] <= 3);
    sum += parts[i];
  }
  assert(sum == 3);
  assert(bic_rank(3, 41, 3, parts, ctx) == 11480);
  bic_ctx_destroy(ctx);

  ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(!bic_params_supported(20, 10, 2, BIC_CACHE_NONE));
  assert(bic_precompute(20, 10, 2, ctx) == 0);
  assert(ctx->comp(20, 10, 2, ctx) == 1);
  bic_run_round_trips(20, 10, 2, 1, NULL, NULL, NULL, NULL, ctx);
  bic_ctx_destroy(ctx);

  for (uint8_t alg = BIC_ALG_PS; alg <= BIC_ALG_AD; ++alg) {
    ctx = bic_ctx_init();
    assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
    assert(bic_ctx_set_unrank_alg(alg, ctx) == 0);
    assert(bic_precompute(20, 10, 2, ctx) != 0);
    bic_ctx_destroy(ctx);
  }
  ctx = bic_ctx_init();
  assert(bic_ctx_set_cache(BIC_CACHE_COMB, ctx) == 0);
  assert(bic_ctx_set_rank_alg(BIC_RANK_ALG_AD, ctx) == 0);
  assert(bic_precompute(20, 10, 2, ctx) != 0);
  bic_ctx_destroy(ctx);

  ctx = bic_ctx_init();
  assert(bic_ctx_set_cache(BIC_CACHE_SMALL_ACC, ctx) == 0);
  assert(bic_params_supported(16, 4, 16, BIC_CACHE_SMALL_ACC));
  assert(!bic_cache_moments_supported(16, 4, 16));
  assert(bic_precompute(16, 4, 16, ctx) != 0);
  bic_ctx_destroy(ctx);

  ctx = bic_ctx_init();
  assert(bic_ctx_set_order(BIC_ORDER_SPIRAL, ctx) == 0);
  assert(bic_precompute(2, 3, 65536, ctx) == 0);
  for (uintx r = 0; r < 6; ++r) {
    uint32_t specialized[3] = {0}, generic[3] = {0};
    assert(bic_ctx_set_unrank_alg(BIC_ALG_DEFAULT, ctx) == 0);
    bic_unrank(specialized, 2, 3, 65536, r, ctx);
    assert(bic_ctx_set_unrank_alg(BIC_ALG_GENERIC, ctx) == 0);
    bic_unrank(generic, 2, 3, 65536, r, ctx);
    for (unsigned i = 0; i < 3; ++i) {
      assert(specialized[i] == generic[i]);
    }
    assert(bic_rank(2, 3, 65536, specialized, ctx) == r);
  }
  bic_ctx_destroy(ctx);

  uint32_t n, d;
  mingen(8, &n, 10, &d);
  assert(n == 4 && d == 2);
  minver(12, &n, 10, &d);
  assert(n == 6 && d <= n);
  assert(bic_params_supported(n, 10, d, BIC_CACHE_NONE));
  assert(compute_bic_no_cache(n, 10, d) >= 4096);
#endif
#if defined(BITINT) && BITINT == 512
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_params_supported(150, 508, 150, BIC_CACHE_NONE));
  const uintx rank = compute_bic_no_cache(150, 507, 150);
  assert(rank < compute_bic_no_cache(150, 508, 150));
  assert(bic_ctx_set_unrank_alg(BIC_ALG_PS, ctx) == 0);
  assert(bic_precompute(150, 508, 150, ctx) == 0);
  uint32_t parts[508] = {0};
  bic_unrank(parts, 150, 508, 150, rank, ctx);
  uint32_t sum = 0;
  for (uint32_t i = 0; i < 508; ++i) {
    assert(parts[i] <= 150);
    sum += parts[i];
  }
  assert(sum == 150);
  assert(bic_rank(150, 508, 150, parts, ctx) == rank);
  bic_ctx_destroy(ctx);
#endif
  return 0;
}
