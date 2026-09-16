#include "api.h"
#include "cache/common.h"
#include "generic/bous.h"
#include "generic/colex.h"
#include "generic/gray.h"
#include "generic/rb_spiral.h"
#include "generic/rbo.h"
#include "generic/spiral.h"
#include "order/bous.h"
#include "order/colex.h"
#include "order/gray.h"
#include "order/rb_spiral.h"
#include "order/rbo.h"
#include "order/spiral.h"
#include "search.h"
#include "strat.h"

typedef struct {
  const char **names;
  const char *target;
} name_search_ctx_t;

const char *bic_rank_alg_names[BIC_RANK_ALG_LENGTH] = {
    "default",
    "generic",
    "ad",
};

uintx compute_bin(const uint32_t n, const uint32_t k, const bic_ctx_t ctx);

uintx compute_bic(const uint32_t n, const uint32_t k, const uint32_t d,
                  const bic_ctx_t ctx);

uint32_t compute_acc(uintx *rop, const uint32_t n, const uint32_t k,
                     const uint32_t d, const bic_ctx_t ctx);

uintx compute_dir(const uint32_t n, const uint32_t k, const uint32_t d,
                  const uint32_t l, const bic_ctx_t ctx);

typedef struct {
  bic_unrank_func_t default_unrank;
  bic_rank_func_t default_rank;
  bic_unrank_func_t generic_unrank;
  bic_rank_func_t generic_rank;
} order_dispatch_entry_t;

static const order_dispatch_entry_t order_dispatch[BIC_ORDER_LENGTH] = {
    {colex_unrank, colex_rank, colex_generic_unrank, colex_generic_rank},
    {gray_unrank, gray_rank, gray_generic_unrank, gray_generic_rank},
    {rbo_unrank, rbo_rank, rbo_generic_unrank, rbo_generic_rank},
    {spiral_unrank, spiral_rank, spiral_generic_unrank, spiral_generic_rank},
    {bous_unrank, bous_rank, bous_generic_unrank, bous_generic_rank},
    {rb_spiral_unrank, rb_spiral_rank, rb_spiral_generic_unrank,
     rb_spiral_generic_rank},
};

static const bic_strategy_func_t strategy_dispatch[BIC_STRATEGY_LENGTH] = {
    mingen,
    minver,
    gen_params_random,
};

static const bic_unrank_func_t colex_unrank_algs[] = {
    NULL,
    colex_unrank_part_sums,
    colex_unrank_acc_linear,
    colex_unrank_acc_bisect,
    colex_unrank_acc_direct,
};

bic_ctx_t bic_ctx_init() {
  bic_ctx_t ctx = (bic_ctx_t)malloc(sizeof(struct bic_ctx_s));
  if (ctx == NULL) {
    return NULL;
  }

  ctx->order = BIC_ORDER_COLEX;
  ctx->unrank_alg = BIC_ALG_DEFAULT;
  ctx->rank_alg = BIC_RANK_ALG_DEFAULT;
  ctx->cache_type = BIC_CACHE_NONE;
  ctx->strategy = BIC_STRATEGY_GEN;
  ctx->scomb_cache_stddev_level = 4;

  ctx->search = strategy_dispatch[BIC_STRATEGY_GEN];
  ctx->unrank = order_dispatch[BIC_ORDER_COLEX].default_unrank;
  ctx->rank = order_dispatch[BIC_ORDER_COLEX].default_rank;

  ctx->bin_cache = NULL;
  ctx->comb_cache = NULL;
  ctx->scomb_cache = NULL;
  ctx->small_acc_cache = NULL;
  ctx->acc_cache = NULL;

  ctx->bin = compute_bin;
  ctx->comp = compute_bic;
  ctx->acc = compute_acc;
  ctx->dir = compute_dir;

  return ctx;
}

void bic_ctx_destroy(bic_ctx_t ctx) { free(ctx); }

bool match_name(const uint32_t i, const void *ctx) {
  const name_search_ctx_t *c = (const name_search_ctx_t *)ctx;
  const uint8_t *s1 = (const uint8_t *)c->names[i];
  const uint8_t *s2 = (const uint8_t *)c->target;

  while (*s1 && (*s1 == *s2)) {
    ++s1;
    ++s2;
  }

  return (*(const uint8_t *)s1 - *(const uint8_t *)s2) == 0;
}

bic_error_t bic_ctx_set_generic_by_name(const char *name,
                                        const find_ctx_t find_ctx,
                                        bic_ctx_t ctx) {
  if (ctx == NULL || name == NULL) {
    return BIC_ERROR_ARGUMENT;
  }
  uint32_t idx;
  name_search_ctx_t c = {.names = find_ctx.names, .target = name};
  linear_search(&idx, 0, find_ctx.length, match_name, &c);
  return find_ctx.setter(idx, ctx);
}

bic_error_t bic_ctx_set_order(const uint8_t order, bic_ctx_t ctx) {
  if (ctx == NULL || order >= BIC_ORDER_LENGTH) {
    return BIC_ERROR_ARGUMENT;
  }

  ctx->unrank = order_dispatch[order].default_unrank;
  ctx->rank = order_dispatch[order].default_rank;
  ctx->order = (bic_order_t)order;
  ctx->unrank_alg = BIC_ALG_DEFAULT;
  ctx->rank_alg = BIC_RANK_ALG_DEFAULT;

  return BIC_SUCCESS;
}

bic_error_t bic_ctx_set_order_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_order_find_ctx, ctx);
}

bic_error_t bic_ctx_set_unrank_alg(const uint8_t alg, bic_ctx_t ctx) {
  if (ctx == NULL || ctx->order >= BIC_ORDER_LENGTH) {
    return BIC_ERROR_ARGUMENT;
  }

  if (alg == BIC_ALG_DEFAULT) {
    ctx->unrank = order_dispatch[ctx->order].default_unrank;
    ctx->unrank_alg = BIC_ALG_DEFAULT;
    return BIC_SUCCESS;
  } else if (alg == BIC_ALG_GENERIC) {
    ctx->unrank = order_dispatch[ctx->order].generic_unrank;
    ctx->unrank_alg = BIC_ALG_GENERIC;
    return BIC_SUCCESS;
  } else if (ctx->order == BIC_ORDER_COLEX) {
    if (alg == BIC_ALG_PS || alg == BIC_ALG_AL || alg == BIC_ALG_AB ||
        alg == BIC_ALG_AD) {
      ctx->unrank = colex_unrank_algs[alg];
      ctx->unrank_alg = (bic_unrank_alg_t)alg;
      return BIC_SUCCESS;
    }
  }

  return BIC_ERROR_ARGUMENT;
}

bic_error_t bic_ctx_set_unrank_alg_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_unrank_alg_find_ctx, ctx);
}

bic_error_t bic_ctx_set_rank_alg(const uint8_t alg, bic_ctx_t ctx) {
  if (ctx == NULL || ctx->order >= BIC_ORDER_LENGTH) {
    return BIC_ERROR_ARGUMENT;
  }

  if (alg == BIC_RANK_ALG_DEFAULT) {
    ctx->rank = order_dispatch[ctx->order].default_rank;
    ctx->rank_alg = BIC_RANK_ALG_DEFAULT;
    return BIC_SUCCESS;
  } else if (alg == BIC_RANK_ALG_GENERIC) {
    ctx->rank = order_dispatch[ctx->order].generic_rank;
    ctx->rank_alg = BIC_RANK_ALG_GENERIC;
    return BIC_SUCCESS;
  } else if (ctx->order == BIC_ORDER_COLEX) {
    if (alg == BIC_RANK_ALG_AD) {
      ctx->rank = colex_rank_acc_direct;
      ctx->rank_alg = BIC_RANK_ALG_AD;
      return BIC_SUCCESS;
    }
  }

  return BIC_ERROR_ARGUMENT;
}

bic_error_t bic_ctx_set_rank_alg_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_rank_alg_find_ctx, ctx);
}

bic_error_t bic_ctx_set_cache(const uint8_t cache, bic_ctx_t ctx) {
  if (ctx == NULL || cache >= BIC_CACHE_LENGTH) {
    return BIC_ERROR_ARGUMENT;
  }

  ctx->cache_type = (bic_cache_t)cache;
  return BIC_SUCCESS;
}

bic_error_t bic_ctx_set_cache_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_cache_find_ctx, ctx);
}

bic_error_t bic_ctx_set_strategy(const uint8_t strategy, bic_ctx_t ctx) {
  if (ctx == NULL || strategy >= BIC_STRATEGY_LENGTH) {
    return BIC_ERROR_ARGUMENT;
  }

  ctx->search = strategy_dispatch[strategy];
  ctx->strategy = (bic_strategy_t)strategy;
  return BIC_SUCCESS;
}

bic_error_t bic_ctx_set_strategy_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_strategy_find_ctx, ctx);
}

bic_error_t bic_ctx_set_defaults(bic_ctx_t ctx) {
  if (ctx == NULL) {
    return BIC_ERROR_ARGUMENT;
  }
  bic_ctx_set_order(BIC_ORDER_COLEX, ctx);
  ctx->cache_type = BIC_CACHE_COMB;
  return BIC_SUCCESS;
}

void bic_run_strategy(const uint32_t m, uint32_t *n, const uint32_t k,
                      uint32_t *d, const bic_ctx_t ctx) {
  ctx->search(m, n, k, d);
}

void bic_unrank(uint32_t *rop, const uint32_t n, const uint32_t k,
                const uint32_t d, const uintx r, const bic_ctx_t ctx) {
  if (k == 0) {
    return;
  }
  ctx->unrank(rop, n, k, d, r, ctx);
}

uintx bic_rank(const uint32_t n, const uint32_t k, const uint32_t d,
               const uint32_t *comp, const bic_ctx_t ctx) {
  if (k == 0) {
    return 0;
  }
  return ctx->rank(n, k, d, comp, ctx);
}

bic_error_t bic_precompute(const uint32_t n, const uint32_t k, const uint32_t d,
                           bic_ctx_t ctx) {
  if (ctx == NULL) {
    return BIC_ERROR_ARGUMENT;
  }
  if (k == 0) {
    return BIC_SUCCESS;
  }
  return build_cache_rec(ctx->cache_type, n, k, d, ctx) ? BIC_ERROR_ALLOCATION
                                                        : BIC_SUCCESS;
}

void bic_free_precomputed(bic_ctx_t ctx) { free_all_caches(ctx); }
