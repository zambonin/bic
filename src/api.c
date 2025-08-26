#include <string.h>

#include "api.h"
#include "cache.h"
#include "colex.h"
#include "gray.h"
#include "rbo.h"
#include "strat.h"
#include "utils.h"

typedef struct {
  const char **names;
  const char *target;
} name_search_ctx_t;

uintx compute_bin(const bic_ctx_t *ctx, const uint32_t n, const uint32_t k);

uintx compute_bic(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                  const uint16_t d);

void compute_acc(uintx *rop, const bic_ctx_t *ctx, const uint16_t n,
                 const uint16_t k, const uint16_t d);

uintx compute_dir(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                  const uint16_t d, const uint16_t l);

bic_ctx_t *bic_ctx_init() {
  bic_ctx_t *ctx = (bic_ctx_t *)malloc(sizeof(bic_ctx_t));
  if (ctx == NULL) {
    return NULL;
  }

  ctx->order = BIC_ORDER_COLEX;
  ctx->unrank_alg = BIC_ALG_DEFAULT;
  ctx->cache_type = BIC_CACHE_NONE;
  ctx->strategy = BIC_STRATEGY_GEN;

  ctx->search = mingen;
  ctx->unrank = colex_unrank;
  ctx->rank = colex_rank;

  ctx->bin_cache = NULL;
  ctx->comb_cache = NULL;
  ctx->scomb_cache = NULL;
  ctx->acc_cache = NULL;

  ctx->bin = compute_bin;
  ctx->comp = compute_bic;
  ctx->acc = compute_acc;
  ctx->dir = compute_dir;

  return ctx;
}

void bic_ctx_destroy(bic_ctx_t *ctx) { free(ctx); }

bool match_name(const uint16_t i, const void *ctx) {
  const name_search_ctx_t *c = (const name_search_ctx_t *)ctx;
  return strcmp(c->names[i], c->target) == 0;
}

uint8_t bic_ctx_set_generic_by_name(bic_ctx_t *ctx, const char *name,
                                    const find_ctx_t find_ctx) {
  uint16_t idx;
  name_search_ctx_t c = {.names = find_ctx.names, .target = name};
  linear_search(&idx, 0, find_ctx.length, match_name, &c);
  return find_ctx.setter(ctx, idx);
}

uint8_t bic_ctx_set_order(bic_ctx_t *ctx, uint8_t order) {
  switch (order) {
  case BIC_ORDER_COLEX:
    ctx->unrank = colex_unrank;
    ctx->rank = colex_rank;
    break;
  case BIC_ORDER_GRAY:
    ctx->unrank = gray_unrank;
    ctx->rank = gray_rank;
    break;
  case BIC_ORDER_RBO:
    ctx->unrank = rbo_unrank;
    ctx->rank = rbo_rank;
    break;
  default:
    return 1;
  }

  ctx->order = (bic_order_t)order;
  return 0;
}

uint8_t bic_ctx_set_order_by_name(bic_ctx_t *ctx, const char *name) {
  return bic_ctx_set_generic_by_name(ctx, name, bic_order_find_ctx);
}

uint8_t bic_ctx_set_unrank_alg(bic_ctx_t *ctx, uint8_t alg) {
  if (ctx->order != BIC_ORDER_COLEX && alg != BIC_ALG_DEFAULT) {
    return 1;
  }

  switch (alg) {
  case BIC_ALG_DEFAULT:
    break;
  case BIC_ALG_PS:
    ctx->unrank = colex_unrank_part_sums;
    break;
  case BIC_ALG_AL:
    ctx->unrank = colex_unrank_acc_linear;
    break;
  case BIC_ALG_AB:
    ctx->unrank = colex_unrank_acc_bisect;
    break;
  case BIC_ALG_AD:
    ctx->unrank = colex_unrank_acc_direct;
    break;
  default:
    return 1;
  }

  ctx->unrank_alg = (bic_unrank_alg_t)alg;
  return 0;
}

uint8_t bic_ctx_set_unrank_alg_by_name(bic_ctx_t *ctx, const char *name) {
  return bic_ctx_set_generic_by_name(ctx, name, bic_unrank_alg_find_ctx);
}

uint8_t bic_ctx_set_cache(bic_ctx_t *ctx, uint8_t cache) {
  if (cache >= BIC_CACHE_LENGTH) {
    return 1;
  }

  ctx->cache_type = (bic_cache_t)cache;
  return 0;
}

uint8_t bic_ctx_set_cache_by_name(bic_ctx_t *ctx, const char *name) {
  return bic_ctx_set_generic_by_name(ctx, name, bic_cache_find_ctx);
}

uint8_t bic_ctx_set_strategy(bic_ctx_t *ctx, uint8_t strategy) {
  switch (strategy) {
  case BIC_STRATEGY_GEN:
    ctx->search = mingen;
    break;
  case BIC_STRATEGY_VER:
    ctx->search = minver;
    break;
  case BIC_STRATEGY_RANDOM:
    ctx->search = gen_params_random;
    break;
  default:
    return 1;
  }

  ctx->strategy = (bic_strategy_t)strategy;
  return 0;
}

uint8_t bic_ctx_set_strategy_by_name(bic_ctx_t *ctx, const char *name) {
  return bic_ctx_set_generic_by_name(ctx, name, bic_strategy_find_ctx);
}

uint8_t bic_ctx_set_defaults(bic_ctx_t *ctx) {
  ctx->order = BIC_ORDER_COLEX;
  ctx->unrank_alg = BIC_ALG_DEFAULT;
  ctx->cache_type = BIC_CACHE_COMB;
  return 0;
}

void bic_run_strategy(bic_ctx_t *ctx, const uint16_t m, uint16_t *n,
                      const uint16_t k, uint16_t *d) {
  ctx->search(m, n, k, d);
}

void bic_unrank(const bic_ctx_t *ctx, uint32_t *rop, const uint16_t n,
                const uint16_t k, const uint16_t d, const uintx r) {
  ctx->unrank(ctx, rop, n, k, d, r);
}

uintx bic_rank(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
               const uint16_t d, const uint32_t *comp) {
  return ctx->rank(ctx, n, k, d, comp);
}

uint8_t bic_precompute(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                       const uint16_t d) {
  return build_cache_rec(ctx, ctx->cache_type, n, k, d);
}

void bic_free_precomputed(bic_ctx_t *ctx) { free_all_caches(ctx); }
