#ifndef API_H
#define API_H

#include "types.h"

typedef uint8_t (*setter_func_t)(bic_ctx_t *, uint8_t);

typedef struct {
  const char **names;
  const setter_func_t setter;
  const uint8_t length;
} find_ctx_t;

bic_ctx_t *bic_ctx_init();

void bic_ctx_destroy(bic_ctx_t *ctx);

uint8_t bic_ctx_set_order(bic_ctx_t *ctx, uint8_t order);

uint8_t bic_ctx_set_order_by_name(bic_ctx_t *ctx, const char *name);

static const char *bic_order_names[BIC_ORDER_LENGTH] = {"colex", "gray", "rbo"};
static const find_ctx_t bic_order_find_ctx = {
    .names = bic_order_names,
    .setter = bic_ctx_set_order,
    .length = BIC_ORDER_LENGTH,
};

uint8_t bic_ctx_set_unrank_alg(bic_ctx_t *ctx, uint8_t alg);

uint8_t bic_ctx_set_unrank_alg_by_name(bic_ctx_t *ctx, const char *name);

static const char *bic_unrank_alg_names[BIC_ALG_LENGTH] = {
    "default", "ps", "al", "ab", "ad",
};
static const find_ctx_t bic_unrank_alg_find_ctx = {
    .names = bic_unrank_alg_names,
    .setter = bic_ctx_set_unrank_alg,
    .length = BIC_ALG_LENGTH,
};

uint8_t bic_ctx_set_cache(bic_ctx_t *ctx, uint8_t cache);

uint8_t bic_ctx_set_cache_by_name(bic_ctx_t *ctx, const char *name);

static const char *bic_cache_names[BIC_CACHE_LENGTH] = {
    "none", "bin", "comb", "scomb", "acc",
};
static const find_ctx_t bic_cache_find_ctx = {
    .names = bic_cache_names,
    .setter = bic_ctx_set_cache,
    .length = BIC_CACHE_LENGTH,
};

uint8_t bic_ctx_set_strategy(bic_ctx_t *ctx, uint8_t strategy);

uint8_t bic_ctx_set_strategy_by_name(bic_ctx_t *ctx, const char *name);

static const char *bic_strategy_names[BIC_STRATEGY_LENGTH] = {
    "gen",
    "ver",
    "random",
};
static const find_ctx_t bic_strategy_find_ctx = {
    .names = bic_strategy_names,
    .setter = bic_ctx_set_strategy,
    .length = BIC_STRATEGY_LENGTH,
};

uint8_t bic_ctx_set_defaults(bic_ctx_t *ctx);

void bic_run_strategy(bic_ctx_t *ctx, const uint16_t m, uint16_t *n,
                      const uint16_t k, uint16_t *d);

void bic_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                const uint16_t k, const uint16_t d, const uintx r);

uintx bic_rank(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
               const uint16_t d, const uint32_t *comp);

uint8_t bic_precompute(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                       const uint16_t d);

void bic_free_precomputed(bic_ctx_t *ctx);

#endif // API_H
