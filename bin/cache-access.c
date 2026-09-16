#include <stdio.h>
#include <string.h>

#include "cache/common.h"
#include "extra.h"

typedef struct {
  uint32_t *access;
  bic_math_bin_func_t bin;
  bic_math_comp_func_t comp;
  bic_math_acc_func_t acc;
  bic_math_dir_func_t dir;
  uint64_t hits;
  uint32_t on_target[BIC_CACHE_LENGTH];
  size_t size;
} access_ctx_t;

static access_ctx_t g_acc_ctx;

uint32_t *fetch(size_t i) { return &(g_acc_ctx.access[i]); }

void count(size_t i) { (*fetch(i))++; }

uint32_t show(size_t i) { return *fetch(i); }

void hit() { (g_acc_ctx.hits)++; }

void on_target(uint8_t cache_type) { (g_acc_ctx.on_target[cache_type])++; }

uintx bin_access(const uint32_t row, const uint32_t col, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_BIN) {
    hit();
    if (col > row) {
      on_target(BIC_CACHE_BIN);
    } else {
      uintx *ptr = bin_cache_get_ptr(row, col, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - ctx->bin_cache->data;
        count(index);
        on_target(BIC_CACHE_BIN);
      }
    }
  }
  return g_acc_ctx.bin(row, col, ctx);
}

uintx comp_access(const uint32_t row, const uint32_t col, const uint32_t d,
                  const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_COMB ||
      ctx->cache_type == BIC_CACHE_SMALL_COMB) {
    hit();
    if (ctx->cache_type == BIC_CACHE_COMB) {
      uintx *ptr = comb_cache_get_ptr(row, col, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - ctx->comb_cache->data;
        count(index);
        on_target(BIC_CACHE_COMB);
      }
    } else if (ctx->cache_type == BIC_CACHE_SMALL_COMB) {
      uintx *ptr = scomb_cache_get_ptr(row, col, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - ctx->scomb_cache->data;
        count(index);
        on_target(BIC_CACHE_SMALL_COMB);
      }
    }
  }
  return g_acc_ctx.comp(row, col, d, ctx);
}

uintx *small_acc_get_val_ptr(const uint32_t j, const uint32_t l,
                             const uint32_t s, const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  if (c == NULL || j >= c->rows || l >= c->cols) {
    return NULL;
  }

  const size_t meta_idx = j * c->cols + l;
  const bic_meta_t m = c->meta[meta_idx];

  if (s < m.base || s >= m.base + m.length) {
    return NULL;
  }

  return &c->data[m.offset + (s - m.base)];
}

uint32_t acc_access(uintx *rop, const uint32_t row, const uint32_t col,
                    const uint32_t d, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC) {
    uint32_t len = 0;
    uintx *ptr = acc_cache_get_ptr(row, col, &len, ctx);
    if (ptr != NULL) {
      const size_t index = ptr - ctx->acc_cache->data;
      for (uint32_t i = 0; i < len; i++) {
        hit();
        count(index + i);
        on_target(BIC_CACHE_ACC);
      }
    }
  }
  return g_acc_ctx.acc(rop, row, col, d, ctx);
}

uintx dir_access(const uint32_t row, const uint32_t col, const uint32_t d,
                 const uint32_t l, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC ||
      ctx->cache_type == BIC_CACHE_SMALL_ACC) {
    if (ctx->cache_type == BIC_CACHE_ACC) {
      if (l == 0) {
        return g_acc_ctx.dir(row, col, d, l, ctx);
      }
      hit();
      uint32_t len = 0;
      uintx *ptr = acc_cache_get_ptr(row, col, &len, ctx);
      if (ptr != NULL && l <= len) {
        const size_t index = ptr - ctx->acc_cache->data;
        count(index + l - 1);
        on_target(BIC_CACHE_ACC);
      }
    } else if (ctx->cache_type == BIC_CACHE_SMALL_ACC) {
      hit();
      uintx *ptr = small_acc_get_val_ptr(col, l, row, ctx);
      if (ptr != NULL) {
        count(ptr - ctx->small_acc_cache->data);
        on_target(BIC_CACHE_SMALL_ACC);
      }
    }
  }
  return g_acc_ctx.dir(row, col, d, l, ctx);
}

static const cli_option_def_t opt_quiet = {
    .opt = {"quiet", required_argument, 0, 'q'},
    .description = "Suppress the cache heatmap when nonzero.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

typedef enum {
  BIC_ACCESS_RANK,
  BIC_ACCESS_UNRANK,
  BIC_ACCESS_ROUND_TRIP,
  BIC_ACCESS_LENGTH,
} bic_access_t;

static const char *access_names[BIC_ACCESS_LENGTH] = {
    "rank",
    "unrank",
    "round-trip",
};

static const char *access_infos[BIC_ACCESS_LENGTH] = {
    "ranking accesses",
    "unranking accesses",
    "both operations (default)",
};

static bic_error_t parse_access(void *value, const char *arg) {
  if (value == NULL || arg == NULL) {
    return BIC_ERROR_ARGUMENT;
  }
  for (uint8_t i = 0; i < BIC_ACCESS_LENGTH; ++i) {
    const int32_t cmp = strcmp(access_names[i], arg);
    if (cmp == 0) {
      *(bic_access_t *)value = (bic_access_t)i;
      return BIC_SUCCESS;
    }
  }
  return BIC_ERROR_ARGUMENT;
}

static const cli_option_def_t opt_operation = {
    .opt = {"operation", required_argument, 0, 'e'},
    .description = "Count rank, unrank, or round-trip accesses (default).",
    .type_description = "operation",
    .parser = parse_access,
    .sub_option_names = access_names,
    .sub_option_infos = access_infos,
    .num_sub_options = BIC_ACCESS_LENGTH,
};

void print_cache_heatmap_bin(const bic_ctx_t ctx) {
  const cache_t *c = ctx->bin_cache;
  for (uint32_t row = 0; row < c->rows; row++) {
    for (uint32_t col = 0; col < c->cols; col++) {
      uintx *ptr = bin_cache_get_ptr(row, col, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - c->data;
        printf("%5d ", show(index));
      } else {
        printf("%5d ", -1);
      }
    }
    printf("\n");
  }
}

void print_cache_heatmap_comb(const bic_ctx_t ctx) {
  const cache_t *c = ctx->comb_cache;
  for (uint32_t row = 0; row < c->rows; row++) {
    for (uint32_t col = 0; col < c->cols; col++) {
      uintx *ptr = comb_cache_get_ptr(row, col, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - c->data;
        printf("%5d ", show(index));
      } else {
        printf("%5d ", -1);
      }
    }
    printf("\n");
  }
}

void print_cache_heatmap_scomb(const bic_ctx_t ctx) {
  const cache_t *c = ctx->scomb_cache;
  for (uint32_t row = 0; row < c->rows; row++) {
    for (uint32_t col = 0; col < c->cols; col++) {
      uintx *ptr = scomb_cache_get_ptr(row, col + 1, ctx);
      if (ptr != NULL) {
        const size_t index = ptr - c->data;
        printf("%5d ", show(index));
      } else {
        printf("%5d ", -1);
      }
    }
    printf("\n");
  }
}

void print_cache_heatmap_acc(const bic_ctx_t ctx) {
  const cache_t *c = ctx->acc_cache;
  const uint32_t columns = c->col_map ? c->col_map_size : c->cols;
  for (uint32_t row = 0; row < c->rows; row++) {
    for (uint32_t col = 0; col < columns; col++) {
      uint32_t len = 0;
      uintx *ptr = acc_cache_get_ptr(row, col, &len, ctx);
      if (ptr == NULL) {
        continue;
      }
      const size_t index = ptr - c->data;

      for (uint32_t depth = 0; depth < c->depth; depth++) {
        if (depth < len) {
          const uint32_t value = show(index + depth);
          printf("%d %d %d %u\n", row, col, depth, value);
        } else {
          printf("%d %d %d %u\n", row, col, depth, 0);
        }
      }
    }
    printf("\n");
  }
}

void print_cache_heatmap_small_acc(const bic_ctx_t ctx) {
  const cache_t *c = ctx->small_acc_cache;
  for (uint32_t j = 1; j < c->rows; j++) {
    for (uint32_t l = 0; l < c->cols; l++) {
      uint32_t lower_s, upper_s;
      if (small_acc_get_bounds(j, l, &lower_s, &upper_s, ctx) != 0) {
        continue;
      }
      for (uint32_t s = lower_s; s < upper_s; ++s) {
        uint32_t val =
            show((size_t)(small_acc_get_val_ptr(j, l, s, ctx) - c->data));
        printf("%u %u %u %u\n", s, j, l, val);
      }
    }
  }
}

void setup_wrappers(bic_ctx_t ctx) {
  g_acc_ctx.bin = ctx->bin;
  ctx->bin = bin_access;

  g_acc_ctx.comp = ctx->comp;
  ctx->comp = comp_access;

  g_acc_ctx.acc = ctx->acc;
  ctx->acc = acc_access;

  g_acc_ctx.dir = ctx->dir;
  ctx->dir = dir_access;
}

static void teardown_wrappers(bic_ctx_t ctx) {
  ctx->bin = g_acc_ctx.bin;
  ctx->comp = g_acc_ctx.comp;
  ctx->acc = g_acc_ctx.acc;
  ctx->dir = g_acc_ctx.dir;
}

static bic_error_t measure_accesses(uint32_t n, uint32_t k, uint32_t d,
                                    uint32_t iterations,
                                    const bic_access_t operation,
                                    bic_ctx_t ctx) {
  uint32_t *parts = (uint32_t *)malloc((size_t)k * sizeof(*parts));
  if (parts == NULL) {
    return BIC_ERROR_ALLOCATION;
  }

  const uintx total = ctx->comp(n, k, d, ctx);
  const bool unrank = operation != BIC_ACCESS_RANK;
  const bool rank = operation != BIC_ACCESS_UNRANK;
  for (uint32_t it = 0; it < iterations; ++it) {
    const uintx r = random_rank_count(total);
    if (unrank) {
      setup_wrappers(ctx);
    }
    bic_unrank(parts, n, k, d, r, ctx);
    if (unrank) {
      teardown_wrappers(ctx);
    }
    if (rank) {
      setup_wrappers(ctx);
    }
    const uintx rr = bic_rank(n, k, d, parts, ctx);
    if (rank) {
      teardown_wrappers(ctx);
    }
    if (r != rr) {
      free(parts);
      return BIC_ERROR_VALIDATION;
    }
    check_valid_bounded_composition(parts, n, k, d);
  }
  free(parts);
  return BIC_SUCCESS;
}

bic_error_t allocate_counters(bic_ctx_t ctx) {
  cache_t *c = NULL;

  switch (ctx->cache_type) {
  case BIC_CACHE_BIN:
    c = ctx->bin_cache;
    break;
  case BIC_CACHE_COMB:
    c = ctx->comb_cache;
    break;
  case BIC_CACHE_SMALL_COMB:
    c = ctx->scomb_cache;
    break;
  case BIC_CACHE_ACC:
    c = ctx->acc_cache;
    break;
  case BIC_CACHE_SMALL_ACC:
    c = ctx->small_acc_cache;
    break;
  default:
    break;
  }

  size_t size = (c != NULL) ? c->length : 0;
  g_acc_ctx.size = size;
  g_acc_ctx.access = (uint32_t *)calloc(size == 0 ? 1 : size, sizeof(uint32_t));
  return g_acc_ctx.access == NULL ? BIC_ERROR_ALLOCATION : BIC_SUCCESS;
}

void print_cache_heatmap(bic_ctx_t ctx) {
  switch (ctx->cache_type) {
  case BIC_CACHE_BIN:
    print_cache_heatmap_bin(ctx);
    break;
  case BIC_CACHE_COMB:
    print_cache_heatmap_comb(ctx);
    break;
  case BIC_CACHE_SMALL_COMB:
    print_cache_heatmap_scomb(ctx);
    break;
  case BIC_CACHE_ACC:
    print_cache_heatmap_acc(ctx);
    break;
  case BIC_CACHE_SMALL_ACC:
    print_cache_heatmap_small_acc(ctx);
    break;
  default:
    break;
  }
}

void free_counters() { free(g_acc_ctx.access); }

int32_t main(int32_t argc, char **argv) {
  bic_error_t error = BIC_SUCCESS;
  uint32_t n = 0;
  uint32_t k = 0;
  uint32_t d = 0;
  uint32_t m = 0;
  char *order = (char *)"colex";
  char *unrank_alg = (char *)"default";
  char *rank_alg = (char *)"default";
  char *cache = (char *)"";
  char *strategy = (char *)"gen";
  bic_access_t operation = BIC_ACCESS_ROUND_TRIP;
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);
  uint32_t stddev = 4;
  uint32_t quiet = 0;

  cli_option_t options[] = {
      {&opt_sum, &n},
      {&opt_parts, &k},
      {&opt_bound, &d},
      {&opt_target, &m},
      {&opt_order, &order},
      {&opt_unrank_alg, &unrank_alg},
      {&opt_rank_alg, &rank_alg},
      {&opt_cache, &cache},
      {&opt_strategy, &strategy},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {&opt_stddev, &stddev},
      {&opt_quiet, &quiet},
      {&opt_operation, &operation},
      {NULL, NULL},
  };

  error = parse_args(argc, argv, options);
  if (error != BIC_SUCCESS) {
    return error;
  }

  srandom(seed);

  long double rate = 0.0;
  uint8_t stddev_level = 0;
  bic_ctx_t ctx = bic_ctx_init();

  if (ctx == NULL) {
    return BIC_ERROR_ALLOCATION;
  }

  if (k == 0 || iterations == 0 || stddev > UINT8_MAX ||
      (uint64_t)n > (uint64_t)k * d || bic_ctx_set_order_by_name(order, ctx) ||
      bic_ctx_set_unrank_alg_by_name(unrank_alg, ctx) ||
      bic_ctx_set_rank_alg_by_name(rank_alg, ctx) ||
      bic_ctx_set_cache_by_name(cache, ctx) ||
      bic_ctx_set_strategy_by_name(strategy, ctx)) {
    print_help(argv[0], options);
    error = BIC_ERROR_ARGUMENT;
    goto cleanup_ctx;
  }

  ctx->scomb_cache_stddev_level = stddev;

  if (m != 0) {
    bic_run_strategy(m, &n, k, &d, ctx);
    if (n == 0 || n == UINT32_MAX) {
      error = BIC_ERROR_RANGE;
      goto cleanup_ctx;
    }
  }

  error = bic_precompute(n, k, d, ctx);
  if (error != BIC_SUCCESS) {
    goto cleanup_precomp;
  }

  error = allocate_counters(ctx);
  if (error != BIC_SUCCESS) {
    goto cleanup_precomp;
  }

  error = measure_accesses(n, k, d, iterations, operation, ctx);
  if (error != BIC_SUCCESS) {
    free_counters();
    goto cleanup_precomp;
  }

  if (g_acc_ctx.hits > 0) {
    rate = (double)g_acc_ctx.on_target[ctx->cache_type] / g_acc_ctx.hits;
  }
  if (ctx->cache_type == BIC_CACHE_SMALL_COMB ||
      ctx->cache_type == BIC_CACHE_SMALL_ACC) {
    stddev_level = ctx->scomb_cache_stddev_level;
  }
  printf("operation = %s, lookups = %llu, cache_hits = %u, "
         "hit_rate = %.6Lf, stddev = %u, size = %zu\n",
         access_names[operation], (unsigned long long)g_acc_ctx.hits,
         g_acc_ctx.on_target[ctx->cache_type], rate, stddev_level,
         g_acc_ctx.size);

  if (!quiet) {
    print_cache_heatmap(ctx);
  }

  free_counters();

cleanup_precomp:
  bic_free_precomputed(ctx);

cleanup_ctx:
  bic_ctx_destroy(ctx);

  return error;
}
