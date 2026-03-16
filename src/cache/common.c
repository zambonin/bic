#include "cache/common.h"
#include "math.h"

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_NONE,
    BIC_CACHE_BIN,  BIC_CACHE_COMB, BIC_CACHE_COMB,
};

uint8_t init_sparse_structures(const uint32_t k, const bic_ctx_t ctx,
                               uint32_t **H_out, size_t *h_size_out,
                               cache_t *c) {
  if (ctx->order != BIC_ORDER_RBO && ctx->order != BIC_ORDER_BOUS) {
    return 0;
  }

  uint16_t L = (uint16_t)lg(k);
  if (L == 0 && k > 1) {
    L = 1;
  }

  uint32_t *H = (uint32_t *)malloc(2 * (size_t)L * sizeof(uint32_t) + 16);
  if (!H) {
    return 1;
  }

  calc_hset(k, L, H, h_size_out);

  if (!cache_allocation_fits((size_t)k + 1, sizeof(int16_t))) {
    free(H);
    return 1;
  }
  c->col_map = (int16_t *)malloc(((size_t)k + 1) * sizeof(int16_t));
  if (!c->col_map) {
    free(H);
    return 1;
  }
  c->col_map_size = k + 1;

  for (uint32_t i = 0; i <= k; ++i) {
    c->col_map[i] = -1;
  }

  *H_out = H;
  return 0;
}

static bool already_built(const bic_cache_t type, const bic_ctx_t ctx) {
  switch (type) {
  case BIC_CACHE_BIN:
    return ctx->bin_cache != NULL;
  case BIC_CACHE_COMB:
    return ctx->comb_cache != NULL;
  case BIC_CACHE_SMALL_COMB:
    return ctx->scomb_cache != NULL;
  case BIC_CACHE_ACC:
    return ctx->acc_cache != NULL;
  case BIC_CACHE_SMALL_ACC:
    return ctx->small_acc_cache != NULL;
  default:
    return true;
  }
}

uint8_t build_cache_rec(const bic_cache_t type, const uint32_t n,
                        const uint32_t k, const uint32_t d, bic_ctx_t ctx) {
  if (type == BIC_CACHE_NONE || already_built(type, ctx)) {
    return 0;
  }

  if (((type == BIC_CACHE_ACC || type == BIC_CACHE_SMALL_ACC ||
        (type == BIC_CACHE_COMB &&
         (ctx->order == BIC_ORDER_RBO || ctx->order == BIC_ORDER_BOUS))) &&
       build_cache_rec(BIC_CACHE_BIN, n, k, d, ctx)) ||
      build_cache_rec(dependencies[type], n, k, d, ctx)) {
    return 1;
  }

  switch (type) {
  case BIC_CACHE_BIN:
    return bin_build_cache(n, k, ctx);
  case BIC_CACHE_COMB:
    return comb_build_cache(n, k, d, ctx);
  case BIC_CACHE_SMALL_COMB:
    return scomb_build_cache(n, k, d, ctx);
  case BIC_CACHE_ACC:
    return acc_build_cache(n, k, d, ctx);
  case BIC_CACHE_SMALL_ACC:
    return small_acc_build_cache(n, k, d, ctx);
  default:
    return 1;
  }

  return 0;
}

void free_all_caches(bic_ctx_t ctx) {
  bin_free_cache(ctx);
  comb_free_cache(ctx);
  scomb_free_cache(ctx);
  acc_free_cache(ctx);
  small_acc_free_cache(ctx);
}
