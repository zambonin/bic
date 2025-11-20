#include <stdio.h>

#include "cache.h"
#include "extra.h"

typedef struct {
  uint32_t *access;
  bic_math_bin_func_t bin;
  bic_math_comp_func_t comp;
  bic_math_acc_func_t acc;
  bic_math_dir_func_t dir;
} access_ctx_t;

static access_ctx_t g_acc_ctx;

uint32_t *fetch(size_t i) { return &(g_acc_ctx.access[i]); }

void count(size_t i) { (*fetch(i))++; }

uint32_t show(size_t i) { return *fetch(i); }

uintx bin_access(const uint32_t row, const uint32_t col, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_BIN) {
    uintx *ptr = bin_cache_get_ptr(row, col, ctx);
    if (ptr != NULL) {
      const size_t index = ptr - ctx->bin_cache->data;
      count(index);
    }
  }
  return g_acc_ctx.bin(row, col, ctx);
}

uintx comp_access(const uint16_t row, const uint16_t col, const uint16_t d,
                  const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_COMB) {
    uintx *ptr = comb_cache_get_ptr(row, col, ctx);
    if (ptr != NULL) {
      const size_t index = ptr - ctx->comb_cache->data;
      count(index);
    }
  } else if (ctx->cache_type == BIC_CACHE_SMALL_COMB) {
    uintx *ptr = scomb_cache_get_ptr(row, col, ctx);
    if (ptr != NULL) {
      const size_t index = ptr - ctx->scomb_cache->data;
      count(index);
    }
  }
  return g_acc_ctx.comp(row, col, d, ctx);
}

uint16_t acc_access(uintx *rop, const uint16_t row, const uint16_t col,
                    const uint16_t d, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC) {
    uint16_t len = 0;
    uintx *ptr = acc_cache_get_ptr(row, col, &len, ctx);
    if (ptr != NULL) {
      const size_t index = ptr - ctx->acc_cache->data;
      for (uint32_t i = 0; i < len; i++) {
        count(index + i);
      }
    }
  }
  return g_acc_ctx.acc(rop, row, col, d, ctx);
}

uintx dir_access(const uint16_t row, const uint16_t col, const uint16_t d,
                 const uint16_t l, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC) {
    uint16_t len = 0;
    uintx *ptr = acc_cache_get_ptr(row, col, &len, ctx);
    if (ptr != NULL && l < len) {
      const size_t index = ptr - ctx->acc_cache->data;
      count(index + l);
    }
  }
  return g_acc_ctx.dir(row, col, d, l, ctx);
}

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
  for (uint32_t row = 0; row < c->rows; row++) {
    for (uint32_t col = 0; col < c->cols; col++) {
      uint16_t len = 0;
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

uint8_t allocate_counters(bic_ctx_t ctx) {
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
  default:
    break;
  }

  size_t size = (c != NULL) ? c->length : 0;
  g_acc_ctx.access = (uint32_t *)calloc(size, sizeof(uint32_t));
  return g_acc_ctx.access == NULL;
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
  default:
    break;
  }
}

void free_counters() { free(g_acc_ctx.access); }

int32_t main(int32_t argc, char **argv) {
  bool error = false;
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint16_t m = 0;
  char *order = (char *)"colex";
  char *algorithm = (char *)"default";
  char *cache = (char *)"";
  char *strategy = (char *)"gen";
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_sum, &n},
      {&opt_parts, &k},
      {&opt_bound, &d},
      {&opt_target, &m},
      {&opt_order, &order},
      {&opt_algorithm, &algorithm},
      {&opt_cache, &cache},
      {&opt_strategy, &strategy},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();

  if (k == 0 || ctx == NULL || (uint64_t)n > (uint64_t)k * d ||
      bic_ctx_set_order_by_name(order, ctx) ||
      bic_ctx_set_unrank_alg_by_name(algorithm, ctx) ||
      bic_ctx_set_cache_by_name(cache, ctx) ||
      bic_ctx_set_strategy_by_name(strategy, ctx)) {
    print_help(argv[0], options);
    error = true;
    goto cleanup_ctx;
  }

  if (m != 0) {
    bic_run_strategy(m, &n, k, &d, ctx);
    if (n == 0 || n == UINT16_MAX) {
      goto cleanup_ctx;
    }
  }

  bic_precompute(n, k, d, ctx);

  if (allocate_counters(ctx)) {
    error = true;
    goto cleanup_precomp;
  }

  setup_wrappers(ctx);

  bic_run_round_trips(n, k, d, iterations, NULL, NULL, NULL, NULL, ctx);

  print_cache_heatmap(ctx);

  free_counters();

cleanup_precomp:
  bic_free_precomputed(ctx);

cleanup_ctx:
  bic_ctx_destroy(ctx);

  return error;
}
