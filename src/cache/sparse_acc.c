#include "cache/common.h"

uintx small_acc_dir_get(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t l, const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || k >= c->rows || l >= c->cols) {
    return compute_dir(n, k, d, l, ctx);
  }

  const bic_meta_t m = c->meta[(size_t)k * c->cols + l];

  if (n < m.base || n >= m.base + m.length) {
    return compute_dir(n, k, d, l, ctx);
  }

  return c->data[m.offset + (n - m.base)];
}

size_t small_acc_cache_length(const uint32_t rows, const uint32_t cols,
                              const uint32_t n, const uint32_t k,
                              const uint32_t d, const uint8_t level,
                              bic_meta_t *const meta, const bic_ctx_t ctx) {
  size_t total = 0;

  for (uint32_t row = 1; row < rows; ++row) {
    for (uint32_t col = 0; col < cols; ++col) {
      const double mean = exp_part_sum(n, k, d, row, col, ctx);
      const double stddev = stddev_part_sum(n, k, d, row, col, ctx);
      const uint32_t lower_bound_n = (col > 0) ? (col - 1) : 0;
      const uint32_t left =
          max(max((int32_t)(mean - level * stddev), 0), lower_bound_n);
      const uint32_t right = min((int32_t)(mean + level * stddev), n);
      const uint32_t len = (right >= left) ? (right - left + 1) : 0;

      if (meta != NULL) {
        bic_meta_t *m = &meta[(size_t)row * cols + col];
        m->base = left;
        m->length = len;
        m->offset = total;
      }

      if (total > SIZE_MAX - len) {
        return SIZE_MAX;
      }
      total += len;
    }
  }

  return total;
}

uint8_t small_acc_build_cache(const uint32_t n, const uint32_t k,
                              const uint32_t d, const bic_ctx_t ctx) {
  if (k == UINT32_MAX) {
    return 1;
  }
  cache_t *c = ctx->small_acc_cache = (cache_t *)calloc(1, sizeof(cache_t));
  if (c == NULL) {
    return 1;
  }

  c->rows = k + 1;
  c->cols = d;
  c->depth = 0;

  if (!cache_product_fits(c->rows, c->cols, sizeof(bic_meta_t))) {
    small_acc_free_cache(ctx);
    return 1;
  }
  c->meta = (bic_meta_t *)calloc((size_t)c->rows * c->cols, sizeof(bic_meta_t));
  if (c->meta == NULL && c->cols != 0) {
    small_acc_free_cache(ctx);
    return 1;
  }
  c->length = small_acc_cache_length(
      c->rows, c->cols, n, k, d, ctx->scomb_cache_stddev_level, c->meta, ctx);
  if (c->length == SIZE_MAX ||
      !cache_allocation_fits(c->length, sizeof(uintx))) {
    small_acc_free_cache(ctx);
    return 1;
  }

  c->data = uintx_alloc(c->length == 0 ? 1 : c->length);
  if (c->data == NULL) {
    small_acc_free_cache(ctx);
    return 1;
  }

  for (uint32_t row = 1; row < c->rows; ++row) {
    for (uint32_t col = 0; col < c->cols; ++col) {
      const bic_meta_t m = c->meta[(size_t)row * c->cols + col];
      if (m.length == 0) {
        continue;
      }
      uintx *const slice = &c->data[m.offset];

      uintx base_sum = 0;
      for (uint32_t i = 0; i < col; ++i) {
        base_sum += ctx->comp(m.base - i, row, d, ctx);
      }
      slice[0] = base_sum;

      for (uint32_t idx = 1; idx < m.length; ++idx) {
        uint32_t n_curr = m.base + idx;
        uintx add = ctx->comp(n_curr, row, d, ctx);
        uintx sub = (n_curr >= col) ? ctx->comp(n_curr - col, row, d, ctx) : 0;
        slice[idx] = slice[idx - 1] + add - sub;
      }
    }
  }

  ctx->dir = small_acc_dir_get;

  return 0;
}

void small_acc_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->small_acc_cache) {
    return;
  }

  free(ctx->small_acc_cache->meta);
  uintx_free(ctx->small_acc_cache->data);
  free(ctx->small_acc_cache);
  ctx->small_acc_cache = NULL;
  ctx->dir = compute_dir;
}

uint8_t small_acc_get_bounds(const uint32_t j, const uint32_t l,
                             uint32_t *lower, uint32_t *upper,
                             const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || j >= c->rows || l >= c->cols) {
    return 1;
  }

  const size_t index = (size_t)j * c->cols + l;
  const bic_meta_t m = c->meta[index];

  *lower = m.base;
  *upper = m.base + m.length;

  return 0;
}
