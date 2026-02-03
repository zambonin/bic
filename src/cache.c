#include "cache.h"
#include "math.h"

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_BIN,
    BIC_CACHE_BIN,  BIC_CACHE_COMB, BIC_CACHE_COMB,
};

static inline bool uses_halving_set(bic_order_t order) {
  return order == BIC_ORDER_RBO || order == BIC_ORDER_BOUS;
}

sparse_cols_t sparse_cols_init(uint32_t k, const bic_ctx_t ctx) {
  sparse_cols_t sc = {
      .columns = NULL,
      .count = 0,
      .col_map = NULL,
      .col_map_size = 0,
      .active = false,
  };

  if (!uses_halving_set(ctx->order)) {
    return sc;
  }

  uint16_t L = (uint16_t)lg(k);
  if (L == 0 && k > 1) {
    L = 1;
  }
  sc.columns = (uint16_t *)malloc(2 * (size_t)L * sizeof(uint16_t) + 16);
  if (!sc.columns) {
    return sc;
  }
  calc_hset((uint16_t)k, L, sc.columns, &sc.count);
  sc.col_map = (int16_t *)malloc(((size_t)k + 1) * sizeof(int16_t));
  if (!sc.col_map) {
    free(sc.columns);
    sc.columns = NULL;
    return sc;
  }
  sc.col_map_size = k + 1;

  for (uint32_t i = 0; i <= k; ++i) {
    sc.col_map[i] = -1;
  }
  for (size_t i = 0; i < sc.count; ++i) {
    sc.col_map[sc.columns[i]] = (int16_t)i;
  }

  sc.active = true;
  return sc;
}

void sparse_cols_filter(sparse_cols_t *sc, uint32_t min_col, uint32_t max_col) {
  if (!sc->active || !sc->columns) {
    return;
  }

  for (size_t i = 0; i < sc->col_map_size; ++i) {
    sc->col_map[i] = -1;
  }
  size_t valid_count = 0;
  for (size_t i = 0; i < sc->count; ++i) {
    uint16_t col = sc->columns[i];
    if (col >= min_col && col < max_col) {
      sc->col_map[col] = (int16_t)valid_count;
      sc->columns[valid_count++] = col;
    }
  }
  sc->count = valid_count;
}

void sparse_cols_transfer(sparse_cols_t *sc, cache_t *c) {
  c->col_map = sc->col_map;
  c->col_map_size = sc->col_map_size;

  sc->col_map = NULL;
  sc->col_map_size = 0;
}

void sparse_cols_free(sparse_cols_t *sc) {
  if (sc->columns) {
    free(sc->columns);
    sc->columns = NULL;
  }
  if (sc->col_map) {
    free(sc->col_map);
    sc->col_map = NULL;
  }
  sc->count = 0;
  sc->col_map_size = 0;
  sc->active = false;
}

static inline uint32_t get_mapped_col(const cache_t *c, uint32_t k) {
  if (c->col_map) {
    if (k >= c->col_map_size || c->col_map[k] == -1) {
      return UINT32_MAX;
    }
    return (uint32_t)c->col_map[k];
  }
  return k;
}

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
  const bic_meta_t m = c->meta[k];
  const uint32_t inner_idx = n - m.base;

  if (inner_idx >= m.length) {
    return NULL;
  }

  return &(c->data[m.offset + inner_idx]);
}

uintx bin_cache_get(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  const uintx *const val = bin_cache_get_ptr(n, k, ctx);
  return (val != NULL) ? *val : 0;
}

uint8_t bin_build_cache(const uint32_t n, const uint32_t k, bic_ctx_t ctx) {
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

size_t comb_cache_length(const uint32_t rows, const uint32_t cols) {
  return rows * cols;
}

uintx *comb_cache_get_ptr(const uint32_t n, const uint32_t k,
                          const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;

  if (n >= c->rows) {
    return NULL;
  }

  uint32_t col_idx = get_mapped_col(c, k);
  if (col_idx == UINT32_MAX) {
    return NULL;
  }

  uint32_t max_col = c->cols;
  if (!c->col_map) {
    if (col_idx >= max_col)
      return NULL;
  }
  if (col_idx >= c->cols)
    return NULL;

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
  cache_t *c = ctx->comb_cache = (cache_t *)malloc(sizeof(cache_t));
  c->col_map = NULL;
  c->col_map_size = 0;

  sparse_cols_t sc = sparse_cols_init(k, ctx);
  uint32_t effective_cols = sc.active ? (uint32_t)sc.count : (k + 1);

  sparse_cols_transfer(&sc, c);

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = 1;
  c->meta = NULL;
  c->length = comb_cache_length(c->rows, c->cols);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    sparse_cols_free(&sc);
    if (c->col_map)
      free(c->col_map);
    free(c);
    ctx->comb_cache = NULL;
    return 1;
  }

  if (sc.active) {
    for (size_t ci = 0; ci < sc.count; ++ci) {
      uint32_t real_k = sc.columns[ci];
      for (uint32_t row = 0; row < c->rows; ++row) {
        *comb_cache_get_ptr(row, real_k, ctx) =
            compute_bic(row, real_k, d, ctx);
      }
    }
  } else {
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
  }

  sparse_cols_free(&sc);
  ctx->comp = comb_cache_get;
  return 0;
}

size_t scomb_cache_length(const uint32_t rows, const uint32_t cols,
                          const uint32_t k, const uint32_t d,
                          const uint8_t level, bic_meta_t *meta,
                          const uint16_t *k_values, const bic_ctx_t ctx) {
  size_t total = 0;
  const uint32_t n = rows - 1;

  for (uint32_t col = 0; col < cols; ++col) {
    const uint32_t kp = k_values ? k_values[col] : (col + 1);

    uint32_t kp_next;
    if (k_values) {
      kp_next = (col + 1 < cols) ? k_values[col + 1] : k;
    } else {
      kp_next = kp + 1;
    }

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

    total += length;
  }

  return total;
}

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

uint8_t scomb_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                          bic_ctx_t ctx) {
  cache_t *c = ctx->scomb_cache = (cache_t *)malloc(sizeof(cache_t));
  c->col_map = NULL;
  c->col_map_size = 0;

  sparse_cols_t sc = sparse_cols_init(k, ctx);

  sparse_cols_filter(&sc, 1, k);

  uint32_t effective_cols = sc.active ? (uint32_t)sc.count : (k - 1);

  sparse_cols_transfer(&sc, c);

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = 1;
  c->meta = (bic_meta_t *)calloc(c->cols, sizeof(bic_meta_t));

  c->length =
      scomb_cache_length(c->rows, c->cols, k, d, ctx->scomb_cache_stddev_level,
                         c->meta, sc.active ? sc.columns : NULL, ctx);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    sparse_cols_free(&sc);
    if (c->col_map)
      free(c->col_map);
    free(c->meta);
    free(c);
    ctx->scomb_cache = NULL;
    return 1;
  }

  if (sc.active) {
    for (uint32_t ci = 0; ci < sc.count; ++ci) {
      const bic_meta_t m = c->meta[ci];
      uint32_t real_k = sc.columns[ci];
      for (uint32_t row = m.base; row < m.base + m.length; ++row) {
        *scomb_cache_get_ptr(row, real_k, ctx) =
            compute_bic(row, real_k, d, ctx);
      }
    }
  } else {
    for (uint32_t col = 0; col < c->cols; ++col) {
      const bic_meta_t m = c->meta[col];
      for (uint32_t row = m.base; row < m.base + m.length; ++row) {
        *scomb_cache_get_ptr(row, col + 1, ctx) =
            compute_bic(row, col + 1, d, ctx);
      }
    }
  }

  sparse_cols_free(&sc);
  ctx->comp = scomb_cache_get;
  return 0;
}

uintx *acc_cache_get_ptr(const uint32_t n, const uint32_t k, uint32_t *len,
                         const bic_ctx_t ctx) {
  const cache_t *c = ctx->acc_cache;
  uint32_t col_idx = get_mapped_col(c, k);

  if (col_idx == UINT32_MAX) {
    return NULL;
  }

  const bic_meta_t m = c->meta[n * c->cols + col_idx];

  if (len != NULL) {
    *len = m.length;
  }

  return &(c->data[m.offset]);
}

uint16_t acc_cache_get(uintx *rop, const uint32_t n, const uint32_t k,
                       const uint32_t d, const bic_ctx_t ctx) {
  (void)d;
  uint32_t len = 0;
  const uintx *source = acc_cache_get_ptr(n, k, &len, ctx);

  for (uint32_t i = 0; i < len; ++i) {
    rop[i] = source[i];
  }

  return len;
}

uintx dir_cache_get(const uint32_t n, const uint32_t k, const uint32_t d,
                    const uint32_t l, const bic_ctx_t ctx) {
  (void)d;
  const uintx *slice = acc_cache_get_ptr(n, k, NULL, ctx);
  return slice[l];
}

size_t acc_cache_length(const uint32_t rows, const uint32_t cols,
                        const uint32_t d, bic_meta_t *const meta) {
  size_t total = 0;

  for (uint32_t row = 0; row < rows; ++row) {
    for (uint32_t col = 0; col < cols; ++col) {
      const uint32_t len = min(row, d) + 2;

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

uint8_t acc_build_cache(const uint32_t n, const uint32_t k, const uint32_t d,
                        bic_ctx_t ctx) {
  cache_t *c = ctx->acc_cache = (cache_t *)malloc(sizeof(cache_t));
  c->col_map = NULL;
  c->col_map_size = 0;

  sparse_cols_t sc = sparse_cols_init(k, ctx);

  sparse_cols_filter(&sc, 0, k);

  uint32_t effective_cols = sc.active ? (uint32_t)sc.count : k;

  sparse_cols_transfer(&sc, c);

  c->rows = n + 1;
  c->cols = effective_cols;
  c->depth = d;
  c->d = d;

  c->meta = (bic_meta_t *)calloc(c->rows * c->cols, sizeof(bic_meta_t));
  c->length = acc_cache_length(c->rows, c->cols, c->d, c->meta);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    sparse_cols_free(&sc);
    if (c->col_map)
      free(c->col_map);
    free(c->meta);
    free(c);
    ctx->acc_cache = NULL;
    return 1;
  }

  uint32_t len = 0;
  for (uint32_t row = 0; row < c->rows; ++row) {
    for (uint32_t ci = 0; ci < c->cols; ++ci) {
      uint32_t real_k = sc.active ? sc.columns[ci] : ci;

      uintx *const slice = acc_cache_get_ptr(row, real_k, &len, ctx);
      slice[0] = 0;
      for (uint32_t l = 0; l < len - 1; ++l) {
        const uintx prev =
            (row < l) ? 0 : comb_cache_get(row - l, real_k, d, ctx);
        slice[l + 1] = slice[l] + prev;
      }
    }
  }

  sparse_cols_free(&sc);
  ctx->acc = acc_cache_get;
  ctx->dir = dir_cache_get;
  return 0;
}

uintx small_acc_dir_get(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t l, const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || l >= c->cols) {
    return compute_dir(n, k, d, l, ctx);
  }

  uint32_t row_idx = get_mapped_col(c, k);
  if (row_idx == UINT32_MAX || row_idx >= c->rows) {
    return compute_dir(n, k, d, l, ctx);
  }

  const bic_meta_t m = c->meta[row_idx * c->cols + l];

  if (n < m.base || n >= m.base + m.length) {
    return compute_dir(n, k, d, l, ctx);
  }

  return c->data[m.offset + (n - m.base)];
}

size_t small_acc_cache_length(const uint32_t rows, const uint32_t cols,
                              const uint32_t n, const uint32_t k,
                              const uint32_t d, const uint8_t level,
                              bic_meta_t *const meta,
                              const uint16_t *k_values, /* NULL for dense */
                              const bic_ctx_t ctx) {
  size_t total = 0;

  if (k_values) {
    for (uint32_t row = 0; row < rows; ++row) {
      uint32_t real_k = k_values[row];

      for (uint32_t col = 0; col < cols; ++col) {
        const double mean = exp_part_sum(n, k, d, real_k, col, ctx);
        const double stddev = stddev_part_sum(n, k, d, real_k, col, ctx);
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
  } else {
    for (uint32_t row = 1; row < rows; ++row) {
      for (uint32_t col = 0; col < cols; ++col) {
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
  }

  return total;
}

uint8_t small_acc_build_cache(const uint32_t n, const uint32_t k,
                              const uint32_t d, const bic_ctx_t ctx) {
  cache_t *c = ctx->small_acc_cache = (cache_t *)malloc(sizeof(cache_t));
  c->col_map = NULL;
  c->col_map_size = 0;

  sparse_cols_t sc = sparse_cols_init(k, ctx);

  sparse_cols_filter(&sc, 1, k + 1);

  uint32_t effective_rows = sc.active ? (uint32_t)sc.count : (k + 1);

  sparse_cols_transfer(&sc, c);

  c->rows = effective_rows;
  c->cols = d;
  c->depth = 0;

  c->meta = (bic_meta_t *)calloc(c->rows * c->cols, sizeof(bic_meta_t));
  c->length = small_acc_cache_length(c->rows, c->cols, n, k, d,
                                     ctx->scomb_cache_stddev_level, c->meta,
                                     sc.active ? sc.columns : NULL, ctx);

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    sparse_cols_free(&sc);
    if (c->col_map)
      free(c->col_map);
    free(c->meta);
    free(c);
    ctx->small_acc_cache = NULL;
    return 1;
  }

  if (sc.active) {
    for (uint32_t ri = 0; ri < c->rows; ++ri) {
      uint32_t real_k = sc.columns[ri];

      for (uint32_t col = 0; col < c->cols; ++col) {
        const bic_meta_t m = c->meta[ri * c->cols + col];
        if (m.length == 0) {
          continue;
        }
        uintx *const slice = &c->data[m.offset];

        uintx base_sum = 0;
        for (uint32_t i = 0; i < col; ++i) {
          base_sum += ctx->comp(m.base - i, real_k, d, ctx);
        }
        slice[0] = base_sum;

        for (uint32_t idx = 1; idx < m.length; ++idx) {
          uint32_t n_curr = m.base + idx;
          uintx add = ctx->comp(n_curr, real_k, d, ctx);
          uintx sub =
              (n_curr >= col) ? ctx->comp(n_curr - col, real_k, d, ctx) : 0;
          slice[idx] = slice[idx - 1] + add - sub;
        }
      }
    }
  } else {
    for (uint32_t row = 1; row < c->rows; ++row) {
      for (uint32_t col = 0; col < c->cols; ++col) {
        const bic_meta_t m = c->meta[row * c->cols + col];
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
          uintx sub =
              (n_curr >= col) ? ctx->comp(n_curr - col, row, d, ctx) : 0;
          slice[idx] = slice[idx - 1] + add - sub;
        }
      }
    }
  }

  sparse_cols_free(&sc);
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

  if (ctx->comb_cache->col_map) {
    free(ctx->comb_cache->col_map);
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

  if (ctx->scomb_cache->col_map) {
    free(ctx->scomb_cache->col_map);
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

  if (ctx->acc_cache->col_map) {
    free(ctx->acc_cache->col_map);
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

  if (ctx->small_acc_cache->col_map) {
    free(ctx->small_acc_cache->col_map);
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
