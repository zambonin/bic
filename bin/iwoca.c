#include <stdio.h>
#include <stdlib.h>

#include "api.h"
#include "extra.h"

#define ITERATIONS 8

typedef struct {
  uint32_t n;
  uint32_t k;
  uint32_t d;
} params_t;

static const params_t PARAMS[] = {
    {68, 136, 1},
    {100, 67, 3},
    {200, 50, 8},
    {525, 25, 42},
};

static const size_t NUM_PARAMS = sizeof(PARAMS) / sizeof(PARAMS[0]);

typedef struct {
  bic_order_t order;
  bic_unrank_alg_t unrank_alg;
  bic_cache_t cache;
  uint8_t sigma;
} perf_config_t;

static const perf_config_t PERF_CONFIGS[] = {
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_NONE, 0},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_BIN, 0},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_COMB, 0},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_ACC, 0},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 2},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 3},
    {BIC_ORDER_COLEX, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 4},

    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_NONE, 0},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_BIN, 0},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_COMB, 0},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_ACC, 0},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_SMALL_COMB, 2},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_SMALL_COMB, 3},
    {BIC_ORDER_COLEX, BIC_ALG_PS, BIC_CACHE_SMALL_COMB, 4},

    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_NONE, 0},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_BIN, 0},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_COMB, 0},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_ACC, 0},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_SMALL_COMB, 2},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_SMALL_COMB, 3},
    {BIC_ORDER_COLEX, BIC_ALG_AB, BIC_CACHE_SMALL_COMB, 4},

    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_NONE, 0},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_BIN, 0},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_COMB, 0},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_ACC, 0},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 2},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 3},
    {BIC_ORDER_RBO, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 4},

    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_NONE, 0},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_BIN, 0},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_COMB, 0},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_ACC, 0},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 2},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 3},
    {BIC_ORDER_SPIRAL, BIC_ALG_GENERIC, BIC_CACHE_SMALL_COMB, 4},
};

static const size_t NUM_PERF_CONFIGS =
    sizeof(PERF_CONFIGS) / sizeof(PERF_CONFIGS[0]);

typedef struct {
  bic_cache_t cache;
  uint8_t sigma;
} storage_config_t;

static const storage_config_t STORAGE_CONFIGS[] = {
    {BIC_CACHE_BIN, 0},        {BIC_CACHE_COMB, 0},
    {BIC_CACHE_ACC, 0},        {BIC_CACHE_SMALL_COMB, 1},
    {BIC_CACHE_SMALL_COMB, 2}, {BIC_CACHE_SMALL_COMB, 3},
    {BIC_CACHE_SMALL_COMB, 4}, {BIC_CACHE_SMALL_ACC, 1},
    {BIC_CACHE_SMALL_ACC, 2},  {BIC_CACHE_SMALL_ACC, 3},
    {BIC_CACHE_SMALL_ACC, 4},
};

static const size_t NUM_STORAGE_CONFIGS =
    sizeof(STORAGE_CONFIGS) / sizeof(STORAGE_CONFIGS[0]);

static size_t get_cache_length(const bic_ctx_t ctx) {
  switch (ctx->cache_type) {
  case BIC_CACHE_BIN:
    return ctx->bin_cache ? ctx->bin_cache->length : 0;
  case BIC_CACHE_COMB:
    return ctx->comb_cache ? ctx->comb_cache->length : 0;
  case BIC_CACHE_SMALL_COMB:
    return ctx->scomb_cache ? ctx->scomb_cache->length : 0;
  case BIC_CACHE_ACC:
    return ctx->acc_cache ? ctx->acc_cache->length : 0;
  case BIC_CACHE_SMALL_ACC:
    return ctx->small_acc_cache ? ctx->small_acc_cache->length : 0;
  default:
    return 0;
  }
}

int main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  printf("type,n,k,d,order,algorithm,cache,sigma,value\n");

  for (size_t p = 0; p < NUM_PARAMS; ++p) {
    const params_t *param = &PARAMS[p];

    for (size_t c = 0; c < NUM_PERF_CONFIGS; ++c) {
      const perf_config_t *cfg = &PERF_CONFIGS[c];

      bic_ctx_set_defaults(ctx);
      bic_ctx_set_order(cfg->order, ctx);
      bic_ctx_set_unrank_alg(cfg->unrank_alg, ctx);
      bic_ctx_set_cache(cfg->cache, ctx);
      ctx->scomb_cache_stddev_level = cfg->sigma;

      if (cfg->cache != BIC_CACHE_NONE) {
        if (bic_precompute(param->n, param->k, param->d, ctx) != 0) {
          continue;
        }
      }

      long double ucycles = 0;
      bic_run_round_trips(param->n, param->k, param->d, ITERATIONS, NULL,
                          &ucycles, NULL, NULL, ctx);

      printf("perf,%u,%u,%u,%s,%s,%s,%u,%.0Lf\n", param->n, param->k, param->d,
             bic_order_names[cfg->order], bic_unrank_alg_names[cfg->unrank_alg],
             bic_cache_names[cfg->cache], cfg->sigma, ucycles / ITERATIONS);

      bic_free_precomputed(ctx);
    }
  }

  static const bic_order_t STORAGE_ORDERS[] = {BIC_ORDER_COLEX, BIC_ORDER_RBO};
  static const size_t NUM_STORAGE_ORDERS =
      sizeof(STORAGE_ORDERS) / sizeof(STORAGE_ORDERS[0]);

  for (size_t p = 0; p < NUM_PARAMS; ++p) {
    const params_t *param = &PARAMS[p];

    for (size_t o = 0; o < NUM_STORAGE_ORDERS; ++o) {
      bic_order_t order = STORAGE_ORDERS[o];

      for (size_t c = 0; c < NUM_STORAGE_CONFIGS; ++c) {
        const storage_config_t *cfg = &STORAGE_CONFIGS[c];

        bic_ctx_set_defaults(ctx);
        bic_ctx_set_order(order, ctx);
        bic_ctx_set_cache(cfg->cache, ctx);
        ctx->scomb_cache_stddev_level = cfg->sigma;

        if (bic_precompute(param->n, param->k, param->d, ctx) != 0) {
          continue;
        }

        size_t elements = get_cache_length(ctx);

        printf("storage,%u,%u,%u,%s,%s,%u,%zu\n", param->n, param->k, param->d,
               bic_order_names[order], bic_cache_names[cfg->cache], cfg->sigma,
               elements);

        bic_free_precomputed(ctx);
      }
    }
  }

  bic_ctx_destroy(ctx);

  return 0;
}
