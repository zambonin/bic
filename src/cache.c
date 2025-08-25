#include "cache.h"
#include "math.h"

#include <string.h>

static const bic_cache_t dependencies[BIC_CACHE_LENGTH] = {
    BIC_CACHE_NONE, BIC_CACHE_NONE, BIC_CACHE_BIN,
    BIC_CACHE_BIN,  BIC_CACHE_COMB,
};

void *cache_get_element(const cache_t *cache, const uint32_t row,
                        const uint32_t col) {
  uint32_t index = row * cache->cols + col;
  size_t offset = index * cache->elem_size;
  return (char *)cache->data + offset;
}

uint8_t generic_setup_cache(cache_t *cache, const uint32_t rows,
                            const uint32_t cols, const size_t elem_size) {
  cache->rows = rows;
  cache->cols = cols;
  cache->elem_size = elem_size;
  cache->total_size = cache->rows * cache->cols * cache->elem_size;
  cache->data = calloc(cache->rows * cache->cols, cache->elem_size);
  return cache->data == NULL;
}

#define GET_CACHE_BIN(ctx, row, col)                                           \
  (*(uintx *)cache_get_element((ctx)->bin_cache, row, col))

uintx bin_from_cache(const bic_ctx_t *ctx, const uint32_t n, const uint32_t k) {
  return (*(uintx *)cache_get_element(ctx->bin_cache, n, k));
}

uint8_t bin_build_cache(bic_ctx_t *ctx, const uint16_t n, const uint16_t k) {
  ctx->bin_cache = (cache_t *)malloc(sizeof(cache_t));
  cache_t *c = ctx->bin_cache;
  if (generic_setup_cache(c, n + k + 1, k, sizeof(uintx))) {
    return 1;
  }

  GET_CACHE_BIN(ctx, 0, 0) = 1;
  for (uint32_t row = 1; row < c->rows; ++row) {
    GET_CACHE_BIN(ctx, row, 0) = 1;
    for (uint16_t col = 1; col <= min(row, c->cols - 1); ++col) {
      GET_CACHE_BIN(ctx, row, col) =
          (uintx)GET_CACHE_BIN(ctx, row - 1, col - 1) +
          (uintx)GET_CACHE_BIN(ctx, row - 1, col);
    }
  }

  ctx->bin = bin_from_cache;

  return 0;
}

#define GET_CACHE_COMB(ctx, row, col)                                          \
  (*(uintx *)cache_get_element((ctx)->comb_cache, row, col))

uintx bic_from_cache_comb(const bic_ctx_t *ctx, const uint16_t n,
                          const uint16_t k, const uint16_t d) {
  (void)d;
  return (*(uintx *)cache_get_element(ctx->comb_cache, n, k));
}

uint8_t comb_build_cache(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                         const uint16_t d) {
  ctx->comb_cache = (cache_t *)malloc(sizeof(cache_t));
  cache_t *c = ctx->comb_cache;
  if (generic_setup_cache(c, n + 1, k + 1, sizeof(uintx))) {
    return 1;
  }

  GET_CACHE_COMB(ctx, 0, 0) = 1;
  for (uint16_t row = 0; row < c->rows; ++row) {
    for (uint16_t col = 1; col < c->cols; ++col) {
      GET_CACHE_COMB(ctx, row, col) = compute_bic(ctx, row, col, d);
    }
  }

  ctx->comp = bic_from_cache_comb;

  return 0;
}

#define GET_CACHE_SCOMB(ctx, row, col)                                         \
  (*(uintx **)cache_get_element((ctx)->scomb_cache, row, col))

uintx bic_from_cache_scomb(const bic_ctx_t *ctx, const uint16_t n,
                           const uint16_t k, const uint16_t d) {
  uintx *row = (*(uintx **)cache_get_element(ctx->scomb_cache, 0, k - 1));
  uint16_t left = (uint16_t)row[0];
  uint16_t right = (uint16_t)row[1];

  if (n < left || n > right) {
    return compute_bic(ctx, n, k, d);
  }

  return row[n - left + 2];
}

uint8_t scomb_build_cache(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                          const uint16_t d) {
  ctx->scomb_cache = (cache_t *)malloc(sizeof(cache_t));
  cache_t *c = ctx->scomb_cache;
  if (generic_setup_cache(ctx->scomb_cache, 1, k - 1, sizeof(uintx *))) {
    return 1;
  }

  double variance = d * (d + 2) / 12;
  uint8_t level = 4;

  for (uint16_t col = 0; col < c->cols; ++col) {
    uint16_t j = col + 1;
    double mean = j * n / k;
    double stddev = asqrt(j * variance * (k - j) / k);
    uint16_t left = max(mean - level * stddev, 0);
    uint16_t right = min(mean + level * stddev, n);

    uint16_t length = 2 + (right - left + 1);
    uintx *part = (uintx *)calloc(length, sizeof(uintx));
    if (part == NULL) {
      return 1;
    }
    c->total_size += length * sizeof(uintx);

    part[0] = left;
    part[1] = right;
    for (uint16_t i = 2; i < length; ++i) {
      part[i] = compute_bic(ctx, left + i - 2, j, d);
    }

    GET_CACHE_SCOMB(ctx, 0, col) = part;
  }

  ctx->comp = bic_from_cache_scomb;

  return 0;
}

#define GET_CACHE_ACC(ctx, row, col)                                           \
  (*(uintx **)cache_get_element((ctx)->acc_cache, row, col))

void acc_from_cache(uintx *rop, const bic_ctx_t *ctx, const uint16_t n,
                    const uint16_t k, const uint16_t d) {
  const size_t length = d + 3;
  uintx *sums = (*(uintx **)cache_get_element(ctx->acc_cache, n, k));
  memcpy(rop, sums, length * sizeof(uintx));
}

uintx dir_from_cache(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                     const uint16_t d, const uint16_t l) {
  const size_t length = d + 3;
  uintx sums[length];
  acc_from_cache(sums, ctx, n, k, d);
  return sums[l];
}

uint8_t acc_build_cache(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                        const uint16_t d) {
  ctx->acc_cache = (cache_t *)malloc(sizeof(cache_t));
  cache_t *c = ctx->acc_cache;
  generic_setup_cache(c, n + 1, k, sizeof(uintx *));
  const size_t length = d + 3;

  for (uint16_t row = 0; row < c->rows; ++row) {
    for (uint16_t col = 0; col < c->cols; ++col) {
      uintx *subarray = (uintx *)calloc(length, sizeof(uintx));
      if (subarray == NULL) {
        return 1;
      }

      compute_acc(subarray, ctx, row, col, d);
      GET_CACHE_ACC(ctx, row, col) = subarray;
      c->total_size += length * sizeof(uintx);
    }
  }

  ctx->acc = acc_from_cache;
  ctx->dir = dir_from_cache;

  return 0;
}

void bin_free_cache(bic_ctx_t *ctx) {
  if (!ctx || !ctx->bin_cache) {
    return;
  }

  free(ctx->bin_cache->data);
  free(ctx->bin_cache);
  ctx->bin_cache = NULL;

  ctx->bin = compute_bin;
}

void comb_free_cache(bic_ctx_t *ctx) {
  if (!ctx || !ctx->comb_cache) {
    return;
  }

  free(ctx->comb_cache->data);
  free(ctx->comb_cache);
  ctx->comb_cache = NULL;

  ctx->comp = compute_bic;
}

void scomb_free_cache(bic_ctx_t *ctx) {
  if (!ctx || !ctx->scomb_cache) {
    return;
  }

  for (uint16_t j = 0; j < ctx->scomb_cache->cols; ++j) {
    free(GET_CACHE_SCOMB(ctx, 0, j));
  }
  free(ctx->scomb_cache->data);
  free(ctx->scomb_cache);
  ctx->scomb_cache = NULL;

  ctx->comp = compute_bic;
}

void acc_free_cache(bic_ctx_t *ctx) {
  if (!ctx || !ctx->acc_cache) {
    return;
  }

  for (uint16_t i = 0; i < ctx->acc_cache->rows; ++i) {
    for (uint16_t j = 0; j < ctx->acc_cache->cols; ++j) {
      free(GET_CACHE_ACC(ctx, i, j));
    }
  }
  free(ctx->acc_cache->data);
  free(ctx->acc_cache);
  ctx->acc_cache = NULL;

  ctx->acc = compute_acc;
  ctx->dir = compute_dir;
}

static bool already_built(const bic_ctx_t *ctx, const bic_cache_t type) {
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

uint8_t build_cache_rec(bic_ctx_t *ctx, bic_cache_t type, const uint16_t n,
                        const uint16_t k, const uint16_t d) {
  if (type == BIC_CACHE_NONE || already_built(ctx, type)) {
    return 0;
  }

  if (build_cache_rec(ctx, dependencies[type], n, k, d)) {
    return 1;
  }

  switch (type) {
  case BIC_CACHE_BIN:
    return bin_build_cache(ctx, n, k);
  case BIC_CACHE_COMB:
    return comb_build_cache(ctx, n, k, d);
  case BIC_CACHE_SMALL_COMB:
    return scomb_build_cache(ctx, n, k, d);
  case BIC_CACHE_ACC:
    return acc_build_cache(ctx, n, k, d);
  default:
    return 1;
  }

  return 0;
}

void free_all_caches(bic_ctx_t *ctx) {
  bin_free_cache(ctx);
  comb_free_cache(ctx);
  scomb_free_cache(ctx);
  acc_free_cache(ctx);
}
