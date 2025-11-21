#include "cache.h"
#include "math.h"

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_BIN,
    BIC_CACHE_BIN,  BIC_CACHE_COMB, BIC_CACHE_COMB,
};

size_t bin_cache_length(const uint16_t rows, const uint16_t cols,
                        bic_meta_t *const meta) {
  size_t total = 0;

  for (uint16_t col = 0; col < cols; ++col) {
    const uint16_t len = rows - col;
    if (meta != NULL) {
      meta[col].offset = total;
      meta[col].length = len;
      meta[col].base = col;
    }
    total += len;
  }

  return total;
}

uintx *bin_cache_get_ptr(const uint16_t n, const uint16_t k,
                         const bic_ctx_t ctx) {
  const cache_t *c = ctx->bin_cache;
  const bic_meta_t m = c->meta[k];
  const uint16_t inner_idx = n - m.base;

  if (inner_idx >= m.length) {
    return NULL;
  }

  return &(c->data[m.offset + inner_idx]);
}

uintx bin_cache_get(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  const uintx *const val = bin_cache_get_ptr(n, k, ctx);
  return (val != NULL) ? *val : 0;
}

uint8_t bin_build_cache(const uint16_t n, const uint16_t k, bic_ctx_t ctx) {
  cache_t *c = ctx->bin_cache = (cache_t *)malloc(sizeof(cache_t));

  c->cols = k + 1;
  c->rows = n + c->cols;
  c->depth = 1;

  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));
  c->length = bin_cache_length(c->rows, c->cols, c->meta);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t row = 0; row < c->rows; ++row) {
    for (uint16_t col = 0; col < min(row + 1, c->cols); ++col) {
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

size_t comb_cache_length(const uint16_t rows, const uint16_t cols) {
  return rows * cols;
}

uintx *comb_cache_get_ptr(const uint16_t n, const uint16_t k,
                          const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;
  if (n >= c->rows || k >= c->cols) {
    return NULL;
  }
  const size_t index = n * c->cols + k;
  return &(c->data[index]);
}

uintx comb_cache_get(const uint16_t n, const uint16_t k, const uint16_t d,
                     const bic_ctx_t ctx) {
  (void)d;
  return *comb_cache_get_ptr(n, k, ctx);
}

uint8_t comb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                         bic_ctx_t ctx) {
  cache_t *c = ctx->comb_cache = (cache_t *)malloc(sizeof(cache_t));

  c->rows = n + 1;
  c->cols = k + 1;
  c->depth = 1;

  c->meta = NULL;
  c->length = comb_cache_length(c->rows, c->cols);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t col = 0; col < c->cols; ++col) {
    *comb_cache_get_ptr(0, col, ctx) = 1;
  }

  for (uint16_t row = 1; row < c->rows; ++row) {
    for (uint16_t col = 1; col < c->cols; ++col) {
      const uintx up = comb_cache_get(row - 1, col, d, ctx);
      const uintx left = comb_cache_get(row, col - 1, d, ctx);
      const uintx inc_enc =
          (row < d + 1) ? 0 : comb_cache_get(row - d - 1, col - 1, d, ctx);
      const uintx res = up + left - inc_enc;
      *comb_cache_get_ptr(row, col, ctx) = res;
    }
  }

  ctx->comp = comb_cache_get;

  return 0;
}

size_t scomb_cache_length(const uint16_t rows, const uint16_t cols,
                          const uint16_t d, const uint8_t level,
                          bic_meta_t *meta, const bic_ctx_t ctx) {
  size_t total = 0;
  const uint16_t n = rows - 1;
  const uint16_t k = cols + 1;

  for (uint16_t col = 0; col < cols; ++col) {
    const double mean = exp_part_sum(n, k, d, col + 1, 0, ctx);
    const double stddev = stddev_part_sum(n, k, d, col + 1, 0, ctx);
    const uint16_t left = max((int32_t)(mean - level * stddev), 0);
    const uint16_t right = min((int32_t)(mean + level * stddev), n);
    const uint16_t length = (right >= left) ? (right - left + 1) : 0;

    if (meta != NULL) {
      meta[col].base = left;
      meta[col].length = length;
      meta[col].offset = total;
    }

    total += length;
  }

  return total;
}

uintx *scomb_cache_get_ptr(const uint16_t n, const uint16_t k,
                           const bic_ctx_t ctx) {
  const cache_t *c = ctx->scomb_cache;
  const bic_meta_t m = c->meta[k - 1];

  if (n < m.base || n >= m.base + m.length) {
    return NULL;
  }

  const size_t index = m.offset + (n - m.base);
  return &(c->data[index]);
}

uintx scomb_cache_get(const uint16_t n, const uint16_t k, const uint16_t d,
                      const bic_ctx_t ctx) {
  const uintx *const val = scomb_cache_get_ptr(n, k, ctx);
  return (val != NULL) ? *val : compute_bic(n, k, d, ctx);
}

uint8_t scomb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                          bic_ctx_t ctx) {
  cache_t *c = ctx->scomb_cache = (cache_t *)malloc(sizeof(cache_t));

  c->rows = n + 1;
  c->cols = k - 1;
  c->depth = 1;

  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));
  c->length = scomb_cache_length(c->rows, c->cols, d,
                                 ctx->scomb_cache_stddev_level, c->meta, ctx);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t col = 0; col < c->cols; ++col) {
    const bic_meta_t m = c->meta[col];
    for (uint16_t row = m.base; row < m.base + m.length; ++row) {
      *scomb_cache_get_ptr(row, col + 1, ctx) =
          compute_bic(row, col + 1, d, ctx);
    }
  }

  ctx->comp = scomb_cache_get;

  return 0;
}

uintx *acc_cache_get_ptr(const uint16_t n, const uint16_t k, uint16_t *len,
                         const bic_ctx_t ctx) {
  const cache_t *c = ctx->acc_cache;
  const bic_meta_t m = c->meta[n * c->cols + k];

  if (len != NULL) {
    *len = m.length;
  }

  return &(c->data[m.offset]);
}

uint16_t acc_cache_get(uintx *rop, const uint16_t n, const uint16_t k,
                       const uint16_t d, const bic_ctx_t ctx) {
  (void)d;
  uint16_t len = 0;
  const uintx *source = acc_cache_get_ptr(n, k, &len, ctx);

  for (uint16_t i = 0; i < len; ++i) {
    rop[i] = source[i];
  }

  return len;
}

uintx dir_cache_get(const uint16_t n, const uint16_t k, const uint16_t d,
                    const uint16_t l, const bic_ctx_t ctx) {
  (void)d;
  const uintx *slice = acc_cache_get_ptr(n, k, NULL, ctx);
  return slice[l];
}

size_t acc_cache_length(const uint16_t rows, const uint16_t cols,
                        const uint16_t d, bic_meta_t *const meta) {
  size_t total = 0;

  for (uint16_t row = 0; row < rows; ++row) {
    for (uint16_t col = 0; col < cols; ++col) {
      const uint16_t len = min(row, d) + 2;

      if (meta != NULL) {
        const size_t i = row * cols + col;
        meta[i].length = len;
        meta[i].offset = total;
        meta[i].base = 0;
      }

      total += len;
    }
  }

  return total;
}

uint8_t acc_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                        bic_ctx_t ctx) {
  cache_t *c = ctx->acc_cache = (cache_t *)malloc(sizeof(cache_t));

  c->rows = n + 1;
  c->cols = k;
  c->depth = d + 2;
  c->d = d;

  c->meta = (bic_meta_t *)calloc(c->rows * c->cols, sizeof(bic_meta_t));
  c->length = acc_cache_length(c->rows, c->cols, c->d, c->meta);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  uint16_t len = 0;
  for (uint16_t row = 0; row < c->rows; ++row) {
    for (uint16_t col = 0; col < c->cols; ++col) {
      uintx *const slice = acc_cache_get_ptr(row, col, &len, ctx);
      slice[0] = 0;
      for (uint16_t l = 0; l < len - 1; ++l) {
        const uintx prev = (row < l) ? 0 : comb_cache_get(row - l, col, d, ctx);
        slice[l + 1] = slice[l] + prev;
      }
    }
  }

  ctx->acc = acc_cache_get;
  ctx->dir = dir_cache_get;

  return 0;
}

uintx small_acc_dir_get(const uint16_t n, const uint16_t k, const uint16_t d,
                        const uint16_t l, const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || k >= c->rows || l >= c->cols) {
    return compute_dir(n, k, d, l, ctx);
  }

  const bic_meta_t m = c->meta[k * c->cols + l];

  if (n < m.base || n >= m.base + m.length) {
    return compute_dir(n, k, d, l, ctx);
  }

  return c->data[m.offset + (n - m.base)];
}

size_t small_acc_cache_length(const uint16_t rows, const uint16_t cols,
                              const uint16_t n, const uint16_t k,
                              const uint16_t d, const uint8_t level,
                              bic_meta_t *const meta, const bic_ctx_t ctx) {
  size_t total = 0;

  for (uint16_t row = 1; row < rows; ++row) {
    for (uint16_t col = 0; col < cols; ++col) {
      const double mean = exp_part_sum(n, k, d, row, col, ctx);
      const double stddev = stddev_part_sum(n, k, d, row, col, ctx);
      const uint16_t lower_bound_n = (col > 0) ? (col - 1) : 0;
      const uint16_t left =
          max(max((int32_t)(mean - level * stddev), 0), lower_bound_n);
      const uint16_t right = min((int32_t)(mean + level * stddev), n);
      const uint16_t len = (right >= left) ? (right - left + 1) : 0;

      if (meta != NULL) {
        bic_meta_t *m = &meta[row * cols + col];
        m->base = left;
        m->length = len;
        m->offset = total;
      }

      total += len;
    }
  }

  return total;
}

uint8_t small_acc_build_cache(const uint16_t n, const uint16_t k,
                              const uint16_t d, const bic_ctx_t ctx) {
  cache_t *c = ctx->small_acc_cache = (cache_t *)malloc(sizeof(cache_t));

  c->rows = k + 1;
  c->cols = d + 1;
  c->depth = 0;

  c->meta = (bic_meta_t *)calloc(c->rows * c->cols, sizeof(bic_meta_t));
  c->length = small_acc_cache_length(
      c->rows, c->cols, n, k, d, ctx->scomb_cache_stddev_level, c->meta, ctx);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t row = 1; row < c->rows; ++row) {
    for (uint16_t col = 0; col < c->cols; ++col) {
      const bic_meta_t m = c->meta[row * c->cols + col];
      if (m.length == 0) {
        continue;
      }
      uintx *const slice = &c->data[m.offset];

      uintx sum = 0;
      for (uint16_t i = 0; i < m.base + 1; ++i) {
        sum += ctx->comp(m.base - i, row, d, ctx);
      }
      slice[0] = sum;

      for (uint16_t l = 1; l < m.length; ++l) {
        slice[l] = slice[l - 1] + ctx->comp(m.base + l, row, d, ctx);
      }
    }
  }

  ctx->dir = small_acc_dir_get;

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

void comb_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->comb_cache) {
    return;
  }

  uintx_free(ctx->comb_cache->data);
  free(ctx->comb_cache);
  ctx->comb_cache = NULL;

  ctx->comp = compute_bic;
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

void acc_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->acc_cache) {
    return;
  }

  free(ctx->acc_cache->meta);
  uintx_free(ctx->acc_cache->data);
  free(ctx->acc_cache);
  ctx->acc_cache = NULL;

  ctx->acc = compute_acc;
  ctx->dir = compute_dir;
}

void small_acc_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->small_acc_cache) {
    return;
  }

  free(ctx->small_acc_cache->meta);
  uintx_free(ctx->small_acc_cache->data);
  free(ctx->small_acc_cache);
  ctx->small_acc_cache = NULL;
}

uint8_t small_acc_get_bounds(const uint16_t j, const uint16_t l,
                             uint32_t *lower, uint32_t *upper,
                             const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || j >= c->rows || l >= c->cols) {
    return 1;
  }

  const size_t index = j * c->cols + l;
  const bic_meta_t m = c->meta[index];

  *lower = m.base;
  *upper = m.base + m.length;

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

uint8_t build_cache_rec(const bic_cache_t type, const uint16_t n,
                        const uint16_t k, const uint16_t d, bic_ctx_t ctx) {
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
    if (ctx->unrank_alg == BIC_ALG_AD) {
      return 0;
    }
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
