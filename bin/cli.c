#include <stdio.h>

#include "extra.h"

int32_t main(int32_t argc, char **argv) {
  bool error = false;
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint16_t m = 0;
  char *order = (char *)"colex";
  char *algorithm = (char *)"default";
  char *cache = (char *)"none";
  char *strategy = (char *)"gen";
  uint32_t iterations = 8;
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
      {&opt_algorithm, &algorithm},
      {&opt_cache, &cache},
      {&opt_strategy, &strategy},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();

  if (k == 0 || ctx == NULL || (uint64_t)n > (uint64_t)k * d ||
      bic_ctx_set_order_by_name(order, ctx) ||
      bic_ctx_set_unrank_alg_by_name(algorithm, ctx) ||
      bic_ctx_set_cache_by_name(cache, ctx) ||
      bic_ctx_set_strategy_by_name(strategy, ctx)) {
    print_help(argv[0], options);
    error = true;
    goto cleanup;
  }

  if (m != 0) {
    bic_run_strategy(m, &n, k, &d, ctx);
    if (n == 0 || n == UINT16_MAX) {
      error = true;
      goto cleanup;
    }
  }

  bic_precompute(n, k, d, ctx);

  bic_run_round_trips(n, k, d, iterations, &utime, &ucycles, &rtime, &rcycles,
                      ctx);

  bic_print_config(n, k, d, ctx);
  bic_print_perf_stats(iterations, utime, ucycles, rtime, rcycles);
  printf("\n");

  bic_free_precomputed(ctx);

cleanup:
  bic_ctx_destroy(ctx);

  return !error;
}
