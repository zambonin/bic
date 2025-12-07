#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cache.h"
#include "extra.h"

static const uint16_t SMALL = 64;

uint16_t uniform_1_to_64(void) { return (random() % SMALL) + 1; }

void gen_small_params_random(uint16_t *n, uint16_t *k, uint16_t *d) {
  *k = uniform_1_to_64();
  *d = uniform_1_to_64();
  *n = (*k * *d + 1) / 2;
}

static const cli_option_def_t opt_stddev = {
    .opt = {"stddev", required_argument, 0, 'v'},
    .description =
        "Configure the maximum amount of standard deviations for the `scomb` "
        "cache to test.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

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

  parse_args(argc, argv, options);

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  bic_ctx_set_cache(BIC_CACHE_COMB, ctx);

  printf("n\tk\td\tstddev\tbin\tcomb\tacc\n");

  for (uint32_t it = 0; it < iterations; ++it) {
    uint16_t n = 0;
    uint16_t k = 0;
    uint16_t d = 0;

    gen_small_params_random(&n, &k, &d, ctx);
    bic_precompute(n, k, d, ctx);

    uint16_t bin_cols = k + 1;
    uint16_t bin_rows = n + bin_cols;

    size_t ragged_size = bin_cache_length(bin_rows, bin_cols, NULL);
    size_t full_size = (size_t)bin_rows * bin_cols;

    uint16_t comb_rows = n + 1;
    uint16_t comb_cols = k + 1;

    size_t comb_size = comb_cache_length(comb_rows, comb_cols);
    uint16_t scomb_rows = n + 1;
    uint16_t scomb_cols = k - 1;

    uint16_t acc_rows = n + 1;
    uint16_t acc_cols = k;

    size_t acc_size = acc_cache_length(acc_rows, acc_cols, d, NULL);
    uint16_t sacc_rows = k + 1;
    uint16_t sacc_cols = d + 1;

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
