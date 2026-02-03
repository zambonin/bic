#ifndef TYPES_H
#define TYPES_H

#include "common.h"

typedef enum {
  BIC_ORDER_COLEX,
  BIC_ORDER_GRAY,
  BIC_ORDER_RBO,
  BIC_ORDER_RBO_MO,
  BIC_ORDER_SPIRAL,
  BIC_ORDER_BOUS,
  BIC_ORDER_LENGTH,
} bic_order_t;

typedef enum {
  BIC_ALG_DEFAULT,
  BIC_ALG_PS,
  BIC_ALG_AL,
  BIC_ALG_AB,
  BIC_ALG_AD,
  BIC_ALG_GENERIC,
  BIC_UNRANK_ALG_LENGTH,
} bic_unrank_alg_t;

typedef enum {
  BIC_CACHE_NONE,
  BIC_CACHE_BIN,
  BIC_CACHE_COMB,
  BIC_CACHE_SMALL_COMB,
  BIC_CACHE_ACC,
  BIC_CACHE_SMALL_ACC,
  BIC_CACHE_LENGTH,
} bic_cache_t;

typedef enum {
  BIC_STRATEGY_GEN,
  BIC_STRATEGY_VER,
  BIC_STRATEGY_RANDOM,
  BIC_STRATEGY_LENGTH,
} bic_strategy_t;

typedef struct {
  size_t offset;
  uint32_t length;
  uint32_t base;
} bic_meta_t;

typedef struct cache_s {
  uintx *data;
  uint32_t rows;
  uint32_t cols;
  int16_t *col_map;
  uint32_t col_map_size;
  uint32_t d;
  uint32_t depth;
  bic_cache_t prereq;
  size_t length;
  bic_meta_t *meta;
} cache_t;

typedef void (*bic_unrank_func_t)(uint32_t *rop, const uint32_t n,
                                  const uint32_t k, const uint32_t d,
                                  const uintx r, const bic_ctx_t ctx);

typedef uintx (*bic_rank_func_t)(const uint32_t n, const uint32_t k,
                                 const uint32_t d, const uint32_t *comp,
                                 const bic_ctx_t ctx);

typedef void (*bic_strategy_func_t)(const uint32_t m, uint32_t *n,
                                    const uint32_t k, uint32_t *d);

typedef uintx (*bic_math_bin_func_t)(const uint32_t n, const uint32_t k,
                                     const bic_ctx_t ctx);

typedef uintx (*bic_math_comp_func_t)(const uint32_t n, const uint32_t k,
                                      const uint32_t d, const bic_ctx_t ctx);

typedef uint16_t (*bic_math_acc_func_t)(uintx *rop, const uint32_t n,
                                        const uint32_t k, const uint32_t d,
                                        const bic_ctx_t ctx);

typedef uintx (*bic_math_dir_func_t)(const uint32_t n, const uint32_t k,
                                     const uint32_t d, const uint32_t l,
                                     const bic_ctx_t ctx);

typedef enum {
  BIC_RANK_ALG_DEFAULT = 0,
  BIC_RANK_ALG_GENERIC = 1,
  BIC_RANK_ALG_AD = 2,
  BIC_RANK_ALG_LENGTH,
} bic_rank_alg_t;

struct bic_ctx_s {
  bic_order_t order;
  bic_unrank_alg_t unrank_alg;
  bic_rank_alg_t rank_alg;
  bic_cache_t cache_type;
  bic_strategy_t strategy;

  bic_unrank_func_t unrank;
  bic_rank_func_t rank;
  bic_strategy_func_t search;

  cache_t *bin_cache;
  cache_t *comb_cache;
  cache_t *scomb_cache;
  cache_t *small_acc_cache;
  cache_t *acc_cache;

  bic_math_bin_func_t bin;
  bic_math_comp_func_t comp;
  bic_math_acc_func_t acc;
  bic_math_dir_func_t dir;

  uint8_t scomb_cache_stddev_level;
};

#endif
