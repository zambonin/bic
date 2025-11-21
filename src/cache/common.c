#include "cache/common.h"
#include "math.h"

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_BIN,
    BIC_CACHE_BIN,  BIC_CACHE_COMB, BIC_CACHE_COMB,
};

uintx *comb_cache_get_value(const uint32_t n, const uint32_t k,
                            const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;
  const size_t index = n * c->cols + k;
  return &(c->data[index]);
}

uintx bic_from_cache_comb(const uint32_t n, const uint32_t k, const uint32_t d,
                          const bic_ctx_t ctx) {
  (void)d;
  return *comb_cache_get_value(n, k, ctx);
}

uint8_t comb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                         bic_ctx_t ctx) {
  ctx->comb_cache = (cache_t *)malloc(sizeof(cache_t));

  cache_t *c = ctx->comb_cache;
  c->rows = n + 1;
  c->cols = k + 1;
  c->depth = 1;
  c->length = c->rows * c->cols * c->depth;
  c->meta = NULL;

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  *comb_cache_get_value(0, 0, ctx) = 1;
  for (uint32_t row = 0; row < c->rows; ++row) {
    for (uint32_t col = 1; col < c->cols; ++col) {
      *comb_cache_get_value(row, col, ctx) = compute_bic(row, col, d, ctx);
    }
  }

  ctx->comp = bic_from_cache_comb;

  return 0;
}

void comb_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->comb_cache) {
    return;
  }

  uintx_free(ctx->comb_cache->data);
  free(ctx->comb_cache);
  ctx->comb_cache = NULL;

  ctx->comp = compute_bic;
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

  if (build_cache_rec(dependencies[type], n, k, d, ctx)) {
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
