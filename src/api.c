#include "api.h"
#include "cache/common.h"
#include "generic/colex.h"
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

  ctx->search = mingen;
  ctx->unrank = colex_unrank;
  ctx->rank = colex_rank;

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

uint8_t bic_ctx_set_generic_by_name(const char *name, const find_ctx_t find_ctx,
                                    bic_ctx_t ctx) {
  uint32_t idx;
  name_search_ctx_t c = {.names = find_ctx.names, .target = name};
  linear_search(&idx, 0, find_ctx.length, match_name, &c);
  return find_ctx.setter(idx, ctx);
}

uint8_t bic_ctx_set_order(const uint8_t order, bic_ctx_t ctx) {
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
  case BIC_ORDER_SPIRAL:
    ctx->unrank = spiral_unrank;
    ctx->rank = spiral_rank;
    break;
  case BIC_ORDER_BOUS:
    ctx->unrank = bous_unrank;
    ctx->rank = bous_rank;
    break;
  case BIC_ORDER_RB_SPIRAL:
    ctx->unrank = rb_spiral_unrank;
    ctx->rank = rb_spiral_rank;
    break;
  default:
    return 1;
  }

  ctx->order = (bic_order_t)order;
  ctx->unrank_alg = BIC_ALG_DEFAULT;
  ctx->rank_alg = BIC_RANK_ALG_DEFAULT;
  return 0;
}

uint8_t bic_ctx_set_order_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_order_find_ctx, ctx);
}

uint8_t bic_ctx_set_unrank_alg(const uint8_t alg, bic_ctx_t ctx) {
  if (ctx->order >= BIC_ORDER_LENGTH ||
      (ctx->order != BIC_ORDER_COLEX && alg != BIC_ALG_DEFAULT &&
       alg != BIC_ALG_GENERIC)) {
    return 1;
  }

  switch (alg) {
  case BIC_ALG_DEFAULT:
    switch (ctx->order) {
    case BIC_ORDER_COLEX:
      ctx->unrank = colex_unrank;
      break;
    case BIC_ORDER_GRAY:
      ctx->unrank = gray_unrank;
      break;
    case BIC_ORDER_RBO:
      ctx->unrank = rbo_unrank;
      break;
    case BIC_ORDER_SPIRAL:
      ctx->unrank = spiral_unrank;
      break;
    case BIC_ORDER_BOUS:
      ctx->unrank = bous_unrank;
      break;
    case BIC_ORDER_RB_SPIRAL:
      ctx->unrank = rb_spiral_unrank;
      break;
    default:
      return 1;
    }
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
  case BIC_ALG_GENERIC:
    switch (ctx->order) {
    case BIC_ORDER_COLEX:
      ctx->unrank = colex_generic_unrank;
      break;
    default:
      return 1;
    }
    break;
  default:
    return 1;
  }

  ctx->unrank_alg = (bic_unrank_alg_t)alg;
  return 0;
}

uint8_t bic_ctx_set_unrank_alg_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_unrank_alg_find_ctx, ctx);
}

uint8_t bic_ctx_set_rank_alg(const uint8_t alg, bic_ctx_t ctx) {
  if (ctx->order >= BIC_ORDER_LENGTH) {
    return 1;
  }
  if (alg == BIC_RANK_ALG_DEFAULT) {
    switch (ctx->order) {
    case BIC_ORDER_COLEX:
      ctx->rank = colex_rank;
      break;
    case BIC_ORDER_GRAY:
      ctx->rank = gray_rank;
      break;
    case BIC_ORDER_RBO:
      ctx->rank = rbo_rank;
      break;
    case BIC_ORDER_SPIRAL:
      ctx->rank = spiral_rank;
      break;
    case BIC_ORDER_BOUS:
      ctx->rank = bous_rank;
      break;
    case BIC_ORDER_RB_SPIRAL:
      ctx->rank = rb_spiral_rank;
      break;
    default:
      return 1;
    }
  } else if (alg == BIC_RANK_ALG_GENERIC) {
    switch (ctx->order) {
    case BIC_ORDER_COLEX:
      ctx->rank = colex_generic_rank;
      break;
    default:
      return 1;
    }
  } else if (alg == BIC_RANK_ALG_AD && ctx->order == BIC_ORDER_COLEX) {
    ctx->rank = colex_rank_acc_direct;
  } else {
    return 1;
  }
  ctx->rank_alg = (bic_rank_alg_t)alg;
  return 0;
}

uint8_t bic_ctx_set_rank_alg_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_rank_alg_find_ctx, ctx);
}

uint8_t bic_ctx_set_cache(const uint8_t cache, bic_ctx_t ctx) {
  if (cache >= BIC_CACHE_LENGTH) {
    return 1;
  }

  ctx->cache_type = (bic_cache_t)cache;
  return 0;
}

uint8_t bic_ctx_set_cache_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_cache_find_ctx, ctx);
}

uint8_t bic_ctx_set_strategy(const uint8_t strategy, bic_ctx_t ctx) {
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

uint8_t bic_ctx_set_strategy_by_name(const char *name, bic_ctx_t ctx) {
  return bic_ctx_set_generic_by_name(name, bic_strategy_find_ctx, ctx);
}

uint8_t bic_ctx_set_defaults(bic_ctx_t ctx) {
  bic_ctx_set_order(BIC_ORDER_COLEX, ctx);
  ctx->cache_type = BIC_CACHE_COMB;
  return 0;
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

uint8_t bic_precompute(const uint32_t n, const uint32_t k, const uint32_t d,
                       bic_ctx_t ctx) {
  if (k == 0) {
    return 0;
  }
  return build_cache_rec(ctx->cache_type, n, k, d, ctx);
}

void bic_free_precomputed(bic_ctx_t ctx) { free_all_caches(ctx); }
