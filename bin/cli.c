#include <stdio.h>

#include "extra.h"

int32_t main(int32_t argc, char **argv) {
  bic_error_t error = BIC_SUCCESS;
  uint32_t n = 0;
  uint32_t k = 0;
  uint32_t d = 0;
  uint32_t m = 0;
  char *order = (char *)"colex";
  char *unrank_alg = (char *)"default";
  char *rank_alg = (char *)"default";
  char *cache = (char *)"none";
  char *strategy = (char *)"gen";
  uint32_t iterations = 8;
  uint32_t stddev = 4;
  uint32_t seed = (uint32_t)time(NULL);

  long double utime = 0;
  long double ucycles = 0;
  long double rtime = 0;
  long double rcycles = 0;

  cli_option_t options[] = {
      {&opt_sum, &n},
      {&opt_parts, &k},
      {&opt_bound, &d},
      {&opt_target, &m},
      {&opt_order, &order},
      {&opt_unrank_alg, &unrank_alg},
      {&opt_rank_alg, &rank_alg},
      {&opt_cache, &cache},
      {&opt_strategy, &strategy},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {&opt_stddev, &stddev},
      {NULL, NULL},
  };

  error = parse_args(argc, argv, options);
  if (error != BIC_SUCCESS) {
    return error;
  }

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();

  if (ctx == NULL) {
    return BIC_ERROR_ALLOCATION;
  }

  if (k == 0 || iterations == 0 || stddev > UINT8_MAX ||
      (uint64_t)n > (uint64_t)k * d || bic_ctx_set_order_by_name(order, ctx) ||
      bic_ctx_set_unrank_alg_by_name(unrank_alg, ctx) ||
      bic_ctx_set_rank_alg_by_name(rank_alg, ctx) ||
      bic_ctx_set_cache_by_name(cache, ctx) ||
      bic_ctx_set_strategy_by_name(strategy, ctx)) {
    print_help(argv[0], options);
    error = BIC_ERROR_ARGUMENT;
    goto cleanup;
  }

  ctx->scomb_cache_stddev_level = stddev;

  if (m != 0) {
    bic_run_strategy(m, &n, k, &d, ctx);
    if (n == 0 || n == UINT32_MAX) {
      error = BIC_ERROR_RANGE;
      goto cleanup;
    }
  }

  error = bic_precompute(n, k, d, ctx);
  if (error != BIC_SUCCESS) {
    goto cleanup;
  }

  error = bic_run_round_trips(n, k, d, iterations, &utime, &ucycles, &rtime,
                              &rcycles, ctx);
  if (error != BIC_SUCCESS) {
    goto cleanup;
  }

  bic_print_config(n, k, d, ctx);
  bic_print_perf_stats(iterations, utime, ucycles, rtime, rcycles);
  printf("\n");

  bic_free_precomputed(ctx);

cleanup:
  bic_ctx_destroy(ctx);

  return error;
}
