#include "cache/common.h"

uintx *acc_cache_get_ptr(const uint32_t n, const uint32_t k, uint32_t *len,
                         const bic_ctx_t ctx) {
  const cache_t *c = ctx->acc_cache;
  if (len != NULL) {
    *len = 0;
  }
  if (c == NULL || n >= c->rows) {
    return NULL;
  }
  uint32_t col_idx = get_mapped_col(c, k);

  if (col_idx >= c->cols) {
    return NULL;
  }

  const bic_meta_t m = ((bic_meta_t *)c->meta)[(size_t)n * c->cols + col_idx];

  if (len != NULL) {
    *len = m.length;
  }

  return m.length == 0 ? NULL : &(c->data[m.offset]);
}

uint32_t acc_cache_get(uintx *rop, const uint32_t n, const uint32_t k,
                       const uint32_t d, const bic_ctx_t ctx) {
  uint32_t len = 0;
  const uintx *source = acc_cache_get_ptr(n, k, &len, ctx);
  if (source == NULL) {
    return compute_acc(rop, n, k, d, ctx);
  }

  for (uint32_t i = 0; i < len; ++i) {
    rop[i] = source[i];
  }

  return len;
}

uintx dir_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t l, const bic_ctx_t ctx) {
  if (l == 0) {
    return 0;
  }
  uint32_t len = 0;
  const uintx *slice = acc_cache_get_ptr(n, k, &len, ctx);
  return slice != NULL && l <= len ? slice[l - 1]
                                   : compute_dir(n, k, d, l, ctx);
}

size_t acc_cache_length(const uint32_t rows, const uint32_t cols,
                        const uint32_t d, bic_meta_t *const meta) {
  size_t total = 0;

  for (uint32_t row = 0; row < rows; ++row) {
    for (uint32_t col = 0; col < cols; ++col) {
      const uint32_t last = min(row, d);
      if (last == UINT32_MAX) {
        return SIZE_MAX;
      }
      const uint32_t len = last + 1;

      if (meta != NULL) {
        const size_t i = (size_t)row * cols + col;
        meta[i].length = len;
        meta[i].offset = total;
        meta[i].base = 0;
      }

      if (total > SIZE_MAX - len) {
        return SIZE_MAX;
      }
      total += len;
    }
  }

  return total;
}

uint8_t acc_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                        bic_ctx_t ctx) {
  if (n == UINT32_MAX) {
    return 1;
  }
  cache_t *c = ctx->acc_cache = (cache_t *)calloc(1, sizeof(cache_t));
  if (c == NULL) {
    return 1;
  }

  uint32_t effective_cols = k;
  uint32_t *H = NULL;
  size_t h_size = 0;

  if (init_sparse_structures(k, ctx, &H, &h_size, c)) {
    acc_free_cache(ctx);
    return 1;
  }

  if (H) {
    size_t valid_count = 0;
    for (size_t i = 0; i < h_size; ++i) {
      if (H[i] < k) {
        c->col_map[H[i]] = (int16_t)valid_count;
        H[valid_count++] = H[i];
      }
    }
    h_size = valid_count;
    effective_cols = (uint32_t)h_size;
  }

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = min(n, d) + 1;
  c->d = d;

  if (!cache_product_fits(c->rows, c->cols, sizeof(bic_meta_t))) {
    free(H);
    acc_free_cache(ctx);
    return 1;
  }
  c->meta = (bic_meta_t *)calloc((size_t)c->rows * c->cols, sizeof(bic_meta_t));
  if (c->meta == NULL) {
    free(H);
    acc_free_cache(ctx);
    return 1;
  }
  c->length = acc_cache_length(c->rows, c->cols, c->d, (bic_meta_t *)c->meta);
  if (c->length == SIZE_MAX ||
      !cache_allocation_fits(c->length, sizeof(uintx))) {
    free(H);
    acc_free_cache(ctx);
    return 1;
  }

  c->data = uintx_alloc(c->length == 0 ? 1 : c->length);
  if (c->data == NULL) {
    free(H);
    acc_free_cache(ctx);
    return 1;
  }

  uint32_t len = 0;
  for (uint32_t row = 0; row < c->rows; ++row) {
    for (uint32_t col = 0; col < c->cols; ++col) {
      uint32_t real_k = (H) ? H[col] : col;

      uintx *const slice = acc_cache_get_ptr(row, real_k, &len, ctx);
      uintx sum = 0;
      for (uint32_t l = 0; l < len; ++l) {
        const uintx prev =
            (row < l) ? 0 : comb_cache_get(row - l, real_k, d, ctx);
        sum += prev;
        slice[l] = sum;
      }
    }
  }

  free(H);

  ctx->acc = acc_cache_get;
  ctx->dir = dir_cache_get;

  return 0;
}

void acc_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->acc_cache) {
    return;
  }

  free(ctx->acc_cache->col_map);
  free(ctx->acc_cache->meta);
  uintx_free(ctx->acc_cache->data);
  free(ctx->acc_cache);
  ctx->acc_cache = NULL;

  ctx->acc = compute_acc;
  ctx->dir = compute_dir;
}
