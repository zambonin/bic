#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cache/common.h"
#include "extra.h"

typedef struct {
  uint64_t hits;
  uint64_t total;
  bic_math_bin_func_t bin;
  bic_math_comp_func_t comp;
  bic_math_acc_func_t acc;
  bic_math_dir_func_t dir;
} hit_ctx_t;

static hit_ctx_t g_hit;

static uintx bin_hit(const uint32_t n, const uint32_t k, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_BIN) {
    g_hit.total++;
    if (bin_cache_get_ptr(n, k, ctx) != NULL) {
      g_hit.hits++;
    }
  }
  return g_hit.bin(n, k, ctx);
}

static uintx comp_hit(const uint32_t n, const uint32_t k, const uint32_t d,
                      const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_COMB) {
    g_hit.total++;
    if (comb_cache_get_ptr(n, k, ctx) != NULL) {
      g_hit.hits++;
    }
  } else if (ctx->cache_type == BIC_CACHE_SMALL_COMB) {
    g_hit.total++;
    if (scomb_cache_get_ptr(n, k, ctx) != NULL) {
      g_hit.hits++;
    }
  }
  return g_hit.comp(n, k, d, ctx);
}

static uint32_t acc_hit(uintx *rop, const uint32_t n, const uint32_t k,
                        const uint32_t d, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC) {
    g_hit.total++;
    uint32_t len = 0;
    if (acc_cache_get_ptr(n, k, &len, ctx) != NULL) {
      g_hit.hits++;
    }
  }
  return g_hit.acc(rop, n, k, d, ctx);
}

static uintx dir_hit(const uint32_t n, const uint32_t k, const uint32_t d,
                     const uint32_t l, const bic_ctx_t ctx) {
  if (ctx->cache_type == BIC_CACHE_ACC) {
    g_hit.total++;
    uint32_t len = 0;
    uintx *ptr = acc_cache_get_ptr(n, k, &len, ctx);
    if (ptr != NULL && l < len) {
      g_hit.hits++;
    }
  } else if (ctx->cache_type == BIC_CACHE_SMALL_ACC) {
    g_hit.total++;
    uint32_t lower, upper;
    if (small_acc_get_bounds(k, l, &lower, &upper, ctx) == 0 && n >= lower &&
        n < upper) {
      g_hit.hits++;
    }
  }
  return g_hit.dir(n, k, d, l, ctx);
}

static void setup_wrappers(bic_ctx_t ctx) {
  g_hit.hits = g_hit.total = 0;
  g_hit.bin = ctx->bin;
  ctx->bin = bin_hit;
  g_hit.comp = ctx->comp;
  ctx->comp = comp_hit;
  g_hit.acc = ctx->acc;
  ctx->acc = acc_hit;
  g_hit.dir = ctx->dir;
  ctx->dir = dir_hit;
}

static void teardown_wrappers(bic_ctx_t ctx) {
  ctx->bin = g_hit.bin;
  ctx->comp = g_hit.comp;
  ctx->acc = g_hit.acc;
  ctx->dir = g_hit.dir;
}

int32_t main(int32_t argc, char **argv) {
  int error = 0;
  uint32_t n = 0;
  uint32_t k = 0;
  uint32_t d = 0;
  uint32_t m = 0;
  char *order = (char *)"colex";
  char *cache = (char *)"";
  char *strategy = (char *)"gen";
  uint32_t iterations = 64;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_sum, &n},
      {&opt_parts, &k},
      {&opt_bound, &d},
      {&opt_target, &m},
      {&opt_order, &order},
      {&opt_cache, &cache},
      {&opt_strategy, &strategy},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  if (parse_args(argc, argv, options) != 0) {
    return EXIT_FAILURE;
  }
  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }
  if (k == 0 || bic_ctx_set_order_by_name(order, ctx) ||
      bic_ctx_set_cache_by_name(cache, ctx) ||
      bic_ctx_set_strategy_by_name(strategy, ctx)) {
    print_help(argv[0], options);
    error = 1;
    goto cleanup;
  }

  if (m != 0) {
    bic_run_strategy(m, &n, k, &d, ctx);
    if (n == 0 || n == UINT32_MAX) {
      error = 1;
      goto cleanup;
    }
  }

  if ((uint64_t)n > (uint64_t)k * d) {
    print_help(argv[0], options);
    error = 1;
    goto cleanup;
  }

  if (bic_precompute(n, k, d, ctx) != 0) {
    fprintf(stderr, "precompute failed for n=%u k=%u d=%u cache=%s\n", n, k, d,
            cache);
    error = 1;
    goto cleanup;
  }

  setup_wrappers(ctx);
  bic_run_round_trips(n, k, d, iterations, NULL, NULL, NULL, NULL, ctx);
  teardown_wrappers(ctx);

  if (g_hit.total == 0) {
    printf("n=%u k=%u d=%u cache=%s: no cache lookups recorded\n", n, k, d,
           cache);
  } else {
    printf("%.4f\n", 100.0 * (double)g_hit.hits / (double)g_hit.total);
  }

cleanup:
  bic_free_precomputed(ctx);
  bic_ctx_destroy(ctx);
  return error;
}
