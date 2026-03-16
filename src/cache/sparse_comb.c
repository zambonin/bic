#include "cache/common.h"

uintx *scomb_cache_get_ptr(const uint32_t n, const uint32_t k,
                           const bic_ctx_t ctx) {
  const cache_t *c = ctx->scomb_cache;
  uint32_t col_idx = get_mapped_col(c, k);
  if (col_idx == UINT32_MAX) {
    return NULL;
  }

  uint32_t meta_idx = (c->col_map) ? col_idx : (k - 1);

  if (meta_idx >= c->cols) {
    return NULL;
  }

  const bic_meta_t m = c->meta[meta_idx];

  if (n < m.base || n >= m.base + m.length) {
    return NULL;
  }

  const size_t index = m.offset + (n - m.base);
  return &(c->data[index]);
}

uintx scomb_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                      const bic_ctx_t ctx) {
  const uintx *const val = scomb_cache_get_ptr(n, k, ctx);
  return (val != NULL) ? *val : compute_bic(n, k, d, ctx);
}

size_t scomb_cache_length(const uint32_t rows, const uint32_t cols,
                          const uint32_t d, const uint8_t level,
                          bic_meta_t *meta, const bic_ctx_t ctx) {
  size_t total = 0;
  const uint32_t n = rows - 1;
  const uint32_t k = cols + 1;

  for (uint32_t col = 0; col < cols; ++col) {
    const uint32_t kp = col + 1;
    const uint32_t kp_next = kp + 1;

    const double mean_target =
        (kp_next >= k) ? (double)n : exp_part_sum(n, k, d, kp_next, 0, ctx);
    const double std_target =
        (kp_next >= k) ? 0.0 : stddev_part_sum(n, k, d, kp_next, 0, ctx);

    const double mean_curr = exp_part_sum(n, k, d, kp, 0, ctx);
    const double std_curr = stddev_part_sum(n, k, d, kp, 0, ctx);

    int32_t left = (int32_t)(mean_curr - level * std_curr);
    int32_t right = (int32_t)(mean_target + level * std_target);

    left = max(left, 0);
    right = max(min(right, (int32_t)n), left);

    const uint32_t length = (right >= left) ? (right - left + 1) : 0;

    if (meta != NULL) {
      meta[col].base = (uint32_t)left;
      meta[col].length = length;
      meta[col].offset = total;
    }

    if (total > SIZE_MAX - length) {
      return SIZE_MAX;
    }
    total += length;
  }

  return total;
}

uint8_t scomb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                          bic_ctx_t ctx) {
  if (n == UINT32_MAX || k == UINT32_MAX) {
    return 1;
  }
  cache_t *c = ctx->scomb_cache = (cache_t *)calloc(1, sizeof(cache_t));
  if (c == NULL) {
    return 1;
  }

  uint32_t effective_cols = k - 1;
  uint32_t *H = NULL;
  size_t h_size = 0;

  if (init_sparse_structures(k, ctx, &H, &h_size, c)) {
    scomb_free_cache(ctx);
    return 1;
  }

  if (H) {
    size_t valid_count = 0;
    for (size_t i = 0; i < h_size; ++i) {
      if (H[i] >= 1 && H[i] < k) {
        c->col_map[H[i]] = (int16_t)valid_count;
        H[valid_count++] = H[i];
      }
    }
    h_size = valid_count;
    effective_cols = (uint32_t)h_size;
  }

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = 1;

  if (!cache_allocation_fits(c->cols, sizeof(bic_meta_t))) {
    free(H);
    scomb_free_cache(ctx);
    return 1;
  }
  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));
  if (c->meta == NULL && c->cols != 0) {
    free(H);
    scomb_free_cache(ctx);
    return 1;
  }

  if (H) {
    size_t total = 0;
    if (!cache_product_fits(c->rows, c->cols, sizeof(uintx))) {
      free(H);
      scomb_free_cache(ctx);
      return 1;
    }
    for (uint32_t col = 0; col < c->cols; ++col) {
      c->meta[col].base = 0;
      c->meta[col].length = c->rows;
      c->meta[col].offset = total;
      total += c->rows;
    }
    c->length = total;
  } else {
    c->length = scomb_cache_length(c->rows, c->cols, d,
                                   ctx->scomb_cache_stddev_level, c->meta, ctx);
  }

  if (c->length == SIZE_MAX ||
      !cache_allocation_fits(c->length, sizeof(uintx))) {
    free(H);
    scomb_free_cache(ctx);
    return 1;
  }
  c->data = uintx_alloc(c->length == 0 ? 1 : c->length);
  if (c->data == NULL) {
    free(H);
    scomb_free_cache(ctx);
    return 1;
  }

  if (H) {
    for (uint32_t col = 0; col < c->cols; ++col) {
      const bic_meta_t m = c->meta[col];
      uint32_t real_k = H[col];
      for (uint32_t row = m.base; row < m.base + m.length; ++row) {
        *scomb_cache_get_ptr(row, real_k, ctx) =
            compute_bic(row, real_k, d, ctx);
      }
    }
    free(H);
  } else {
    for (uint32_t col = 0; col < c->cols; ++col) {
      const bic_meta_t m = c->meta[col];
      for (uint32_t row = m.base; row < m.base + m.length; ++row) {
        *scomb_cache_get_ptr(row, col + 1, ctx) =
            compute_bic(row, col + 1, d, ctx);
      }
    }
  }

  ctx->comp = scomb_cache_get;

  return 0;
}

void scomb_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->scomb_cache) {
    return;
  }

  free(ctx->scomb_cache->col_map);
  free(ctx->scomb_cache->meta);
  uintx_free(ctx->scomb_cache->data);
  free(ctx->scomb_cache);
  ctx->scomb_cache = NULL;

  ctx->comp = compute_bic;
}
