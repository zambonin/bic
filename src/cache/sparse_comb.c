#include "cache/common.h"

uintx *scomb_cache_get_ptr(const uint32_t n, const uint32_t k,
                           const bic_ctx_t ctx) {
  const cache_t *c = ctx->scomb_cache;
  uint32_t meta_idx = k - 1;

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
    const double mean_curr = exp_part_sum(n, k, d, kp, 0, ctx);
    const double std_curr = stddev_part_sum(n, k, d, kp, 0, ctx);

    int32_t left = (int32_t)(mean_curr - level * std_curr);
    int32_t right = (int32_t)(mean_curr + level * std_curr);

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

  c->rows = n + 1;
  c->cols = k - 1;
  c->depth = 1;

  if (!cache_allocation_fits(c->cols, sizeof(bic_meta_t))) {
    scomb_free_cache(ctx);
    return 1;
  }
  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));
  if (c->meta == NULL && c->cols != 0) {
    scomb_free_cache(ctx);
    return 1;
  }

  c->length = scomb_cache_length(c->rows, c->cols, d,
                                 ctx->scomb_cache_stddev_level, c->meta, ctx);

  if (c->length == SIZE_MAX ||
      !cache_allocation_fits(c->length, sizeof(uintx))) {
    scomb_free_cache(ctx);
    return 1;
  }
  c->data = uintx_alloc(c->length == 0 ? 1 : c->length);
  if (c->data == NULL) {
    scomb_free_cache(ctx);
    return 1;
  }

  for (uint32_t col = 0; col < c->cols; ++col) {
    const bic_meta_t m = c->meta[col];
    for (uint32_t row = m.base; row < m.base + m.length; ++row) {
      *scomb_cache_get_ptr(row, col + 1, ctx) =
          compute_bic(row, col + 1, d, ctx);
    }
  }

  ctx->comp = scomb_cache_get;

  return 0;
}

void scomb_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->scomb_cache) {
    return;
  }

  free(ctx->scomb_cache->meta);
  uintx_free(ctx->scomb_cache->data);
  free(ctx->scomb_cache);
  ctx->scomb_cache = NULL;

  ctx->comp = compute_bic;
}
