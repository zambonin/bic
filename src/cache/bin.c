#include "cache/common.h"
#include "math.h"

size_t bin_cache_length(const uint32_t rows, const uint32_t cols,
                        bic_meta_t *const meta) {
  size_t total = 0;

  for (uint32_t col = 0; col < cols; ++col) {
    const uint32_t len = rows - col;
    if (meta != NULL) {
      meta[col].offset = total;
      meta[col].length = len;
      meta[col].base = col;
    }
    total += len;
  }

  return total;
}

uintx *bin_cache_get_ptr(const uint32_t n, const uint32_t k,
                         const bic_ctx_t ctx) {
  const cache_t *c = ctx->bin_cache;
  if (k >= c->cols || n >= c->rows || k > n) {
    return NULL;
  }
  const bic_meta_t m = ((bic_meta_t *)c->meta)[k];
  const uint32_t inner_idx = n - m.base;

  if (inner_idx >= m.length) {
    return NULL;
  }

  return &(c->data[m.offset + inner_idx]);
}

uintx bin_cache_get(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  if (k > n) {
    return 0;
  }
  const uintx *const val = bin_cache_get_ptr(n, k, ctx);
  return (val != NULL) ? *val : compute_bin(n, k, ctx);
}

uint8_t bin_build_cache(const uint32_t n, const uint32_t k, bic_ctx_t ctx) {
  if ((uint64_t)n + (uint64_t)k + 1 > UINT32_MAX) {
    return 1;
  }
  cache_t *c = ctx->bin_cache = (cache_t *)calloc(1, sizeof(cache_t));
  if (c == NULL) {
    return 1;
  }

  c->cols = k + 1;
  c->rows = n + c->cols;
  c->depth = 1;

  if (!cache_allocation_fits(c->cols, sizeof(bic_meta_t)) ||
      !cache_product_fits(c->rows, c->cols, sizeof(uintx))) {
    bin_free_cache(ctx);
    return 1;
  }
  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));
  if (c->meta == NULL) {
    bin_free_cache(ctx);
    return 1;
  }
  c->length = bin_cache_length(c->rows, c->cols, (bic_meta_t *)c->meta);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    bin_free_cache(ctx);
    return 1;
  }

  for (uint32_t row = 0; row < c->rows; ++row) {
    for (uint32_t col = 0; col < min(row + 1, c->cols); ++col) {
      if (col == 0 || col == row) {
        *bin_cache_get_ptr(row, col, ctx) = 1;
      } else {
        const uintx up = bin_cache_get(row - 1, col - 1, ctx);
        const uintx left = bin_cache_get(row - 1, col, ctx);
        const uintx res = up + left;
        *bin_cache_get_ptr(row, col, ctx) = res;
      }
    }
  }

  ctx->bin = bin_cache_get;

  return 0;
}

void bin_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->bin_cache) {
    return;
  }

  free(ctx->bin_cache->meta);
  uintx_free(ctx->bin_cache->data);
  free(ctx->bin_cache);
  ctx->bin_cache = NULL;

  ctx->bin = compute_bin;
}
