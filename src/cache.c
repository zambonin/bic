#include "cache.h"
#include "math.h"

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_BIN,
    BIC_CACHE_BIN,  BIC_CACHE_COMB,
};

uintx bin_from_cache(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  if (k > n) {
    return 0;
  }

  const cache_t *c = ctx->bin_cache;
  const bin_cache_meta_t *meta = (bin_cache_meta_t *)c->meta;
  const size_t index = meta->offsets[k] + (n - k);

  return c->data[index];
}

uint8_t bin_build_cache(const uint16_t n, const uint16_t k, bic_ctx_t ctx) {
  ctx->bin_cache = (cache_t *)malloc(sizeof(cache_t));
  bin_cache_meta_t *meta = (bin_cache_meta_t *)malloc(sizeof(bin_cache_meta_t));

  cache_t *c = ctx->bin_cache;
  c->rows = 1;
  c->cols = k + 1;
  c->depth = 1;
  c->length = 0;
  c->meta = meta;

  meta->max_row = n + c->cols;
  meta->offsets = (size_t *)calloc(c->cols, sizeof(size_t));
  for (uint16_t col = 0; col < c->cols; ++col) {
    meta->offsets[col] = c->length;
    c->length += (meta->max_row - col);
  }

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t row = 0; row < meta->max_row; ++row) {
    for (uint16_t col = 0; col < min(row + 1, c->cols); ++col) {
      const size_t index = meta->offsets[col] + (row - col);
      if (col == 0 || col == row) {
        c->data[index] = 1;
      } else {
        c->data[index] = bin_from_cache(row - 1, col - 1, ctx) +
                         bin_from_cache(row - 1, col, ctx);
      }
    }
  }

  ctx->bin = bin_from_cache;

  return 0;
}

uintx *comb_cache_get_value(const uint32_t n, const uint32_t k,
                            const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;
  const size_t index = n * c->cols + k;
  return &(c->data[index]);
}

uintx bic_from_cache_comb(const uint16_t n, const uint16_t k, const uint16_t d,
                          const bic_ctx_t ctx) {
  (void)d;
  return *comb_cache_get_value(n, k, ctx);
}

uint8_t comb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
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

  for (uint16_t col = 0; col < c->cols; ++col) {
    *comb_cache_get_value(0, col, ctx) = 1;
  }

  for (uint16_t row = 1; row < c->rows; ++row) {
    for (uint16_t col = 1; col < c->cols; ++col) {
      const uintx up = bic_from_cache_comb(row - 1, col, d, ctx);
      const uintx left = bic_from_cache_comb(row, col - 1, d, ctx);
      const uintx inc_enc =
          (row < d + 1) ? 0 : bic_from_cache_comb(row - d - 1, col - 1, d, ctx);
      const uintx res = up + left - inc_enc;
      *comb_cache_get_value(row, col, ctx) = res;
    }
  }

  ctx->comp = bic_from_cache_comb;

  return 0;
}

uintx bic_from_cache_scomb(const uint16_t n, const uint16_t k, const uint16_t d,
                           const bic_ctx_t ctx) {
  const cache_t *c = ctx->scomb_cache;
  const scomb_cache_meta_t m = ((scomb_cache_meta_t *)c->meta)[k - 1];

  if (n < m.left || n > m.right) {
    return compute_bic(n, k, d, ctx);
  }

  return c->data[m.offset + (n - m.left)];
}

uint8_t scomb_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                          bic_ctx_t ctx) {
  ctx->scomb_cache = (cache_t *)malloc(sizeof(cache_t));

  cache_t *c = ctx->scomb_cache;
  c->rows = n + 1;
  c->cols = k - 1;
  c->depth = 1;
  c->length = 0;

  scomb_cache_meta_t *meta =
      (scomb_cache_meta_t *)calloc(c->cols, sizeof(scomb_cache_meta_t));
  c->meta = meta;

  const uint8_t level = ctx->scomb_cache_stddev_level;
  for (uint16_t col = 0; col < c->cols; ++col) {
    uint16_t j = col + 1;
    double mean = exp_part_sum(n, k, d, j, 0, ctx);
    double stddev = stddev_part_sum(n, k, d, j, 0, ctx);
    meta[col].left = max((int32_t)(mean - level * stddev), 0);
    meta[col].right = min((int32_t)(mean + level * stddev), n);
    meta[col].offset = c->length;
    c->length += (meta[col].right - meta[col].left + 1);
  }

  c->data = uintx_alloc(c->length);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t col = 0; col < c->cols; ++col) {
    scomb_cache_meta_t m = meta[col];
    uint16_t j = col + 1;
    uint16_t length = m.right - m.left + 1;
    for (uint16_t i = 0; i < length; ++i) {
      c->data[m.offset + i] = compute_bic(m.left + i, j, d, ctx);
    }
  }

  ctx->comp = bic_from_cache_scomb;

  return 0;
}

uintx *acc_cache_get_value(const uint16_t n, const uint16_t k, const uint16_t d,
                           const uint16_t l, const bic_ctx_t ctx) {
  (void)d;
  const cache_t *c = ctx->acc_cache;
  const size_t index = (n * c->cols * c->depth) + (k * c->depth) + l;
  return &(c->data[index]);
}

uint16_t acc_from_cache(uintx *rop, const uint16_t n, const uint16_t k,
                        const uint16_t d, const bic_ctx_t ctx) {
  const cache_t *c = ctx->acc_cache;
  const acc_cache_meta_t *meta = (acc_cache_meta_t *)c->meta;
  const uintx *source = acc_cache_get_value(n, k, d, 0, ctx);
  for (uint32_t i = 0; i < c->depth; ++i) {
    rop[i] = source[i];
  }
  return meta->lengths[n * c->cols + k];
}

uintx dir_from_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                     const uint16_t l, const bic_ctx_t ctx) {
  return *acc_cache_get_value(n, k, d, l, ctx);
}

uint8_t acc_build_cache(const uint16_t n, const uint16_t k, const uint16_t d,
                        bic_ctx_t ctx) {
  ctx->acc_cache = (cache_t *)malloc(sizeof(cache_t));
  acc_cache_meta_t *meta = (acc_cache_meta_t *)malloc(sizeof(acc_cache_meta_t));

  cache_t *c = ctx->acc_cache;
  c->rows = n + 1;
  c->cols = k;
  c->depth = d + 2;
  c->meta = meta;

  meta->lengths = (uint16_t *)malloc(c->rows * c->cols * sizeof(uint16_t));
  c->data = uintx_alloc(c->rows * c->cols * c->depth);
  if (c->data == NULL) {
    return 1;
  }

  for (uint16_t row = 0; row < c->rows; ++row) {
    for (uint16_t col = 0; col < c->cols; ++col) {
      uintx *dest = acc_cache_get_value(row, col, d, 0, ctx);
      uint16_t len = compute_acc(dest, row, col, d, ctx);
      meta->lengths[row * c->cols + col] = len;
    }
  }

  ctx->acc = acc_from_cache;
  ctx->dir = dir_from_cache;

  return 0;
}

void bin_free_cache(bic_ctx_t ctx) {
  if (!ctx || !ctx->bin_cache) {
    return;
  }

  bin_cache_meta_t *meta = (bin_cache_meta_t *)ctx->bin_cache->meta;
  free(meta->offsets);
  free(meta);
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

  acc_cache_meta_t *meta = (acc_cache_meta_t *)ctx->acc_cache->meta;
  free(meta->lengths);
  free(meta);
  uintx_free(ctx->acc_cache->data);
  free(ctx->acc_cache);
  ctx->acc_cache = NULL;

  ctx->acc = compute_acc;
  ctx->dir = compute_dir;
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
}
