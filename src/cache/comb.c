#include "cache/common.h"

size_t comb_cache_length(const uint32_t rows, const uint32_t cols) {
  return (size_t)rows * cols;
}

uintx *comb_cache_get_ptr(const uint32_t n, const uint32_t k,
                          const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;

  if (n >= c->rows) {
    return NULL;
  }

  const uint32_t col_idx = get_mapped_col(c, k);
  if (col_idx >= c->cols) {
    return NULL;
  }

  const size_t index = (size_t)n * c->cols + col_idx;
  return &(c->data[index]);
}

uintx comb_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                     const bic_ctx_t ctx) {
  const uintx *ptr = comb_cache_get_ptr(n, k, ctx);
  return (ptr != NULL) ? *ptr : compute_bic(n, k, d, ctx);
}

uint8_t comb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                         bic_ctx_t ctx) {
  if (n == UINT32_MAX || k == UINT32_MAX) {
    return 1;
  }
  cache_t *c = ctx->comb_cache = (cache_t *)calloc(1, sizeof(cache_t));
  if (c == NULL) {
    return 1;
  }

  uint32_t effective_cols = k + 1;
  uint32_t *H = NULL;
  size_t h_size = 0;

  if (init_sparse_structures(k, ctx, &H, &h_size, c)) {
    comb_free_cache(ctx);
    return 1;
  }

  if (H) {
    for (size_t i = 0; i < h_size; ++i) {
      c->col_map[H[i]] = (int16_t)i;
    }
    effective_cols = (uint32_t)h_size;
  }

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = 1;

  c->meta = NULL;
  if (!cache_product_fits(c->rows, c->cols, sizeof(uintx))) {
    free(H);
    comb_free_cache(ctx);
    return 1;
  }
  c->length = comb_cache_length(c->rows, c->cols);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    free(H);
    comb_free_cache(ctx);
    return 1;
  }

  if (H) {
    for (size_t j = 0; j < h_size; ++j) {
      uint32_t real_k = H[j];
      for (uint32_t row = 0; row < c->rows; ++row) {
        *comb_cache_get_ptr(row, real_k, ctx) =
            compute_bic(row, real_k, d, ctx);
      }
    }
    free(H);
  } else {
    for (uint32_t col = 0; col < c->cols; ++col) {
      *comb_cache_get_ptr(0, col, ctx) = 1;
    }

    for (uint32_t row = 1; row < c->rows; ++row) {
      for (uint32_t col = 1; col < c->cols; ++col) {
        const uintx up = comb_cache_get(row - 1, col, d, ctx);
        const uintx left = comb_cache_get(row, col - 1, d, ctx);
        const uintx inc_enc =
            (row < d + 1) ? 0 : comb_cache_get(row - d - 1, col - 1, d, ctx);
        const uintx res = up - inc_enc + left;
        *comb_cache_get_ptr(row, col, ctx) = res;
      }
    }
  }

  ctx->comp = comb_cache_get;

  return 0;
}

void comb_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->comb_cache) {
    return;
  }

  free(ctx->comb_cache->col_map);
  uintx_free(ctx->comb_cache->data);
  free(ctx->comb_cache);
  ctx->comb_cache = NULL;

  ctx->comp = compute_bic;
}
