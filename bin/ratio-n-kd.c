#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "extra.h"

static const uint32_t K = 32;
static const uint32_t D = 255;

int32_t run_ratio_sweep(uint32_t iterations) {
  int32_t error = 1;
  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  bic_ctx_set_cache(BIC_CACHE_COMB, ctx);
  if (bic_precompute(K * D, K, D, ctx) != 0) {
    goto cleanup;
  }

  for (double r = 0.05; r <= 2.001; r += 0.05) {
    uint32_t n = (uint32_t)(r * ((double)K * D / 2.0));
    if (n > K * D) {
      n = K * D;
    }

    uintx *ranks = uintx_alloc(iterations);
    uint32_t *comp = (uint32_t *)calloc(K, sizeof(uint32_t));
    if (ranks == NULL || comp == NULL) {
      free(comp);
      uintx_free(ranks);
      goto cleanup;
    }
    const uintx count = ctx->comp(n, K, D, ctx);
    for (uint32_t i = 0; i < iterations; ++i) {
      ranks[i] = random_rank_count(count);
    }

    for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
      if (bic_ctx_set_order(i, ctx) != 0) {
        free(comp);
        uintx_free(ranks);
        goto cleanup;
      }
      long double total_cycles = 0;

      for (uint32_t it = 0; it < iterations; ++it) {
        uintx r_in = ranks[it];
        long double t = 0, c = 0;
        PERF(t, c, bic_unrank(comp, n, K, D, r_in, ctx), unrank);
        (void)t;
        total_cycles += c;
      }
      printf("n=%u, k=%u, d=%u, ratio=%.2f, o=%s, "
             "unrank_avg_cyc=%.2Lf\n",
             n, K, D, r, bic_order_names[i], total_cycles / iterations);
    }

    free(comp);
    uintx_free(ranks);
  }

  error = 0;

cleanup:
  bic_free_precomputed(ctx);
  bic_ctx_destroy(ctx);
  return error;
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 1000;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  if (parse_args(argc, argv, options) != 0 || iterations == 0) {
    return EXIT_FAILURE;
  }
  srandom(seed);

  return run_ratio_sweep(iterations);
}
