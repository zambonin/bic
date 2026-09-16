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
  bic_error_t error = BIC_SUCCESS;
  uint32_t iterations = 10;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
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
  bic_ctx_set_cache(BIC_CACHE_COMB, ctx);

  printf("n\tk\td\tcomb_dense\tcomb_sparse\tcomb_save%%\tacc_dense\tacc_"
         "sparse\tacc_save%%\tscomb_dense\tscomb_sparse\tscomb_save%%\n");

  for (uint32_t it = 0; it < iterations; ++it) {
    uint32_t n = 0;
    uint32_t k = 0;
    uint32_t d = 0;

    bool supported = false;
    for (uint32_t tries = 0; tries < 1024 && !supported; ++tries) {
      gen_small_params_random(&n, &k, &d);
      supported = bic_cache_moments_supported(n, k, d);
    }
    error = supported ? bic_precompute(n, k, d, ctx) : BIC_ERROR_RANGE;
    if (error != BIC_SUCCESS) {
      break;
    }

    uint32_t H[64];
    size_t sparse_cols = 0;
    calc_hset(k, (uint16_t)lg(k), H, &sparse_cols);

    const size_t comb_dense = comb_cache_length(n + 1, k + 1);
    const size_t comb_sparse = comb_cache_length(n + 1, sparse_cols);
    const size_t acc_dense = acc_cache_length(n + 1, k, d, NULL);
    const size_t acc_sparse = acc_cache_length(n + 1, sparse_cols, d, NULL);
    const size_t scomb_dense = scomb_cache_length(
        n + 1, k - 1, d, ctx->scomb_cache_stddev_level, NULL, ctx);
    const size_t scomb_sparse = comb_cache_length(n + 1, sparse_cols);

    if (comb_dense == SIZE_MAX || comb_sparse == SIZE_MAX ||
        acc_dense == SIZE_MAX || acc_sparse == SIZE_MAX ||
        scomb_dense == SIZE_MAX || scomb_sparse == SIZE_MAX) {
      error = BIC_ERROR_ALLOCATION;
      bic_free_precomputed(ctx);
      break;
    }

    double comb_save =
        comb_dense > 0 ? 100.0 * (1.0 - (double)comb_sparse / comb_dense) : 0.0;
    double acc_save =
        acc_dense > 0 ? 100.0 * (1.0 - (double)acc_sparse / acc_dense) : 0.0;
    double scomb_save = scomb_dense > 0
                            ? 100.0 * (1.0 - (double)scomb_sparse / scomb_dense)
                            : 0.0;

    printf("%d\t%d\t%d\t%zu\t%zu\t%.2f\t%zu\t%zu\t%.2f\t%zu\t%zu\t%.2f\n", n, k,
           d, comb_dense, comb_sparse, comb_save, acc_dense, acc_sparse,
           acc_save, scomb_dense, scomb_sparse, scomb_save);
    bic_free_precomputed(ctx);
  }

  bic_ctx_destroy(ctx);

  return error;
}
