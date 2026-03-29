#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cache/common.h"
#include "extra.h"

static const uint32_t SMALL = 64;

uint32_t uniform_1_to_64(void) { return (random() % SMALL) + 1; }

void gen_small_params_random(uint32_t *n, uint32_t *k, uint32_t *d) {
  *k = uniform_1_to_64();
  *d = uniform_1_to_64();
  *n = (*k * *d + 1) / 2;
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);
  uint32_t max_stddev = 1;

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {&opt_stddev, &max_stddev},
      {NULL, NULL},
  };

  if (parse_args(argc, argv, options) != 0) {
    return EXIT_FAILURE;
  }

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  bic_ctx_set_cache(BIC_CACHE_COMB, ctx);

  printf("n\tk\td\tstddev\tbin\tcomb\tacc\n");

  for (uint32_t it = 0; it < iterations; ++it) {
    uint32_t n = 0;
    uint32_t k = 0;
    uint32_t d = 0;

    gen_small_params_random(&n, &k, &d);
    bic_precompute(n, k, d, ctx);

    uint32_t bin_cols = k + 1;
    uint32_t bin_rows = n + bin_cols;

    size_t ragged_size = bin_cache_length(bin_rows, bin_cols, NULL);
    size_t full_size = (size_t)bin_rows * bin_cols;

    uint32_t comb_rows = n + 1;
    uint32_t comb_cols = k + 1;

    size_t comb_size = comb_cache_length(comb_rows, comb_cols);
    uint32_t scomb_rows = n + 1;
    uint32_t scomb_cols = k - 1;

    uint32_t acc_rows = n + 1;
    uint32_t acc_cols = k;

    size_t acc_size = acc_cache_length(acc_rows, acc_cols, d, NULL);
    uint32_t sacc_rows = k + 1;
    uint32_t sacc_cols = d + 1;

    for (uint32_t s = 1; s <= max_stddev; ++s) {
      ctx->scomb_cache_stddev_level = s;

      size_t scomb_size_val = scomb_cache_length(
          scomb_rows, scomb_cols, d, ctx->scomb_cache_stddev_level, NULL, ctx);
      size_t sacc_size_val =
          small_acc_cache_length(sacc_rows, sacc_cols, n, k, d,
                                 ctx->scomb_cache_stddev_level, NULL, ctx);

      double ragged_save_percent = (1.0 - (double)ragged_size / full_size);
      double scomb_save_percent = (1.0 - (double)scomb_size_val / comb_size);
      double sacc_save_percent = (1.0 - (double)sacc_size_val / acc_size);

      printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", n, k, d, s,
             ragged_save_percent, scomb_save_percent, sacc_save_percent);
    }
    bic_free_precomputed(ctx);
  }

  bic_ctx_destroy(ctx);

  return 0;
}
