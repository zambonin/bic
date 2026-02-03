#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "extra.h"

static const uint32_t K = 32;
static const uint32_t D = 255;

void run_ratio_sweep(uint32_t iterations) {
  bic_ctx_t ctx = bic_ctx_init();
  if (!ctx)
    return;

  bic_ctx_set_cache(BIC_CACHE_COMB, ctx);
  bic_precompute(K * D, K, D, ctx);

  printf("ratio n k d");
  for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
    printf(" %s", bic_order_names[i]);
  }
  printf("\n");

  for (double r = 0.05; r <= 2.001; r += 0.05) {
    uint32_t n = (uint32_t)(r * ((double)K * D / 2.0));
    if (n > K * D)
      n = K * D;

    long double avg_cycles[BIC_ORDER_LENGTH];
    for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
      avg_cycles[i] = 0;
    }

    uintx *ranks = (uintx *)calloc(iterations, sizeof(uintx));
    for (uint32_t i = 0; i < iterations; ++i) {
      ranks[i] = random_rank(n, K, D);
    }
    uint32_t *comp = (uint32_t *)calloc(K, sizeof(uint32_t));

    for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
      bic_ctx_set_order(i, ctx);
      long double total_cycles = 0;

      for (uint32_t it = 0; it < iterations; ++it) {
        uintx r_in = ranks[it];
        long double t = 0, c = 0;
        PERF(t, c, bic_unrank(comp, n, K, D, r_in, ctx), unrank);
        total_cycles += c;
      }
      avg_cycles[i] = total_cycles / iterations;
    }

    printf("%.2f %u %u %u", r, n, K, D);
    for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
      printf(" %.2Lf", avg_cycles[i]);
    }
    printf("\n");

    fprintf(stderr, "Ratio %.2f (n=%u):", r, n);
    for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
      fprintf(stderr, " %s=%.0Lf", bic_order_names[i], avg_cycles[i]);
      if (i < BIC_ORDER_LENGTH - 1) {
        fprintf(stderr, ",");
      }
    }
    fprintf(stderr, "\n");

    free(comp);
    free(ranks);
  }

  bic_free_precomputed(ctx);
  bic_ctx_destroy(ctx);
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 1000;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);
  srandom(seed);

  run_ratio_sweep(iterations);

  return 0;
}
