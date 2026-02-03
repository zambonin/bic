#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cache.h"
#include "extra.h"

static const uint32_t SMALL = 64;

uint32_t uniform_1_to_64(void) { return (random() % SMALL) + 1; }

void gen_small_params_random(uint32_t *n, uint32_t *k, uint32_t *d) {
  *k = uniform_1_to_64();
  *d = uniform_1_to_64();
  *n = (*k * *d + 1) / 2;
}

size_t get_cache_size(bic_ctx_t ctx, bic_cache_t type, uint32_t n, uint32_t k,
                      uint32_t d) {
  if (build_cache_rec(type, n, k, d, ctx))
    return 0;

  size_t size = 0;
  switch (type) {
  case BIC_CACHE_COMB:
    if (ctx->comb_cache)
      size = ctx->comb_cache->length;
    break;
  case BIC_CACHE_ACC:
    if (ctx->acc_cache)
      size = ctx->acc_cache->length;
    break;
  case BIC_CACHE_SMALL_COMB:
    if (ctx->scomb_cache)
      size = ctx->scomb_cache->length;
    break;
  default:
    break;
  }
  free_all_caches(ctx);
  return size;
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 10;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  printf("n\tk\td\tcomb_dense\tcomb_sparse\tcomb_save%%\tacc_dense\tacc_"
         "sparse\tacc_save%%\tscomb_dense\tscomb_sparse\tscomb_save%%\n");

  for (uint32_t it = 0; it < iterations; ++it) {
    uint32_t n = 0;
    uint32_t k = 0;
    uint32_t d = 0;

    gen_small_params_random(&n, &k, &d);

    bic_ctx_set_order(BIC_ORDER_COLEX, ctx);
    size_t comb_dense = get_cache_size(ctx, BIC_CACHE_COMB, n, k, d);
    bic_ctx_set_order(BIC_ORDER_RBO, ctx);
    size_t comb_sparse = get_cache_size(ctx, BIC_CACHE_COMB, n, k, d);

    bic_ctx_set_order(BIC_ORDER_COLEX, ctx);
    size_t acc_dense = get_cache_size(ctx, BIC_CACHE_ACC, n, k, d);
    bic_ctx_set_order(BIC_ORDER_RBO, ctx);
    size_t acc_sparse = get_cache_size(ctx, BIC_CACHE_ACC, n, k, d);

    bic_ctx_set_order(BIC_ORDER_COLEX, ctx);
    size_t scomb_dense = get_cache_size(ctx, BIC_CACHE_SMALL_COMB, n, k, d);
    bic_ctx_set_order(BIC_ORDER_RBO, ctx);
    size_t scomb_sparse = get_cache_size(ctx, BIC_CACHE_SMALL_COMB, n, k, d);

    double comb_save = 100.0 * (1.0 - (double)comb_sparse / comb_dense);
    double acc_save = 100.0 * (1.0 - (double)acc_sparse / acc_dense);
    double scomb_save = 100.0 * (1.0 - (double)scomb_sparse / scomb_dense);

    printf("%d\t%d\t%d\t%zu\t%zu\t%.2f\t%zu\t%zu\t%.2f\t%zu\t%zu\t%.2f\n", n, k,
           d, comb_dense, comb_sparse, comb_save, acc_dense, acc_sparse,
           acc_save, scomb_dense, scomb_sparse, scomb_save);
  }

  bic_ctx_destroy(ctx);

  return 0;
}
