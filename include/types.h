#ifndef TYPES_H
#define TYPES_H

#include "common.h"

typedef enum {
  BIC_ORDER_COLEX,
  BIC_ORDER_GRAY,
  BIC_ORDER_RBO,
  BIC_ORDER_LENGTH,
} bic_order_t;

typedef enum {
  BIC_ALG_DEFAULT,
  BIC_ALG_PS,
  BIC_ALG_AL,
  BIC_ALG_AB,
  BIC_ALG_AD,
  BIC_ALG_LENGTH,
} bic_unrank_alg_t;

typedef enum {
  BIC_CACHE_NONE,
  BIC_CACHE_BIN,
  BIC_CACHE_COMB,
  BIC_CACHE_SMALL_COMB,
  BIC_CACHE_ACC,
  BIC_CACHE_LENGTH,
} bic_cache_t;

typedef enum {
  BIC_STRATEGY_GEN,
  BIC_STRATEGY_VER,
  BIC_STRATEGY_RANDOM,
  BIC_STRATEGY_LENGTH,
} bic_strategy_t;

typedef struct cache_s {
  void *data;
  uint32_t rows;
  uint32_t cols;
  size_t elem_size;
  char *name;
  uint8_t type;
  size_t total_size;
  bic_cache_t prereq;
} cache_t;

typedef void (*bic_unrank_func_t)(const bic_ctx_t *, uint32_t *, const uint16_t,
                                  const uint16_t, const uint16_t, const uintx);

typedef uintx (*bic_rank_func_t)(const bic_ctx_t *, const uint16_t,
                                 const uint16_t, const uint16_t,
                                 const uint32_t *);

typedef void (*bic_strategy_func_t)(const uint16_t, uint16_t *, const uint16_t,
                                    uint16_t *);

typedef uintx (*bic_math_bin_func_t)(const bic_ctx_t *ctx, const uint32_t n,
                                     const uint32_t k);

typedef uintx (*bic_math_comp_func_t)(const bic_ctx_t *ctx, const uint16_t n,
                                      const uint16_t k, const uint16_t d);

typedef void (*bic_math_acc_func_t)(uintx *rop, const bic_ctx_t *ctx,
                                    const uint16_t n, const uint16_t k,
                                    const uint16_t d);

typedef uintx (*bic_math_dir_func_t)(const bic_ctx_t *ctx, const uint16_t n,
                                     const uint16_t k, const uint16_t d,
                                     const uint16_t l);

typedef struct bic_ctx_s {
  bic_order_t order;
  bic_unrank_alg_t unrank_alg;
  bic_cache_t cache_type;
  bic_strategy_t strategy;

  bic_unrank_func_t unrank;
  bic_rank_func_t rank;
  bic_strategy_func_t search;

  cache_t *bin_cache;
  cache_t *comb_cache;
  cache_t *scomb_cache;
  cache_t *acc_cache;

  bic_math_bin_func_t bin;
  bic_math_comp_func_t comp;
  bic_math_acc_func_t acc;
  bic_math_dir_func_t dir;
} bic_ctx_t;

#endif // TYPES_H
