#include <stdio.h>

#include "api.h"
#include "extra.h"

#define IS_BIT_SET(mask, idx) ((mask) & (1 << (idx)))

typedef struct {
  uint8_t order;
  uint8_t unrank_alg;
  uint8_t rank_alg;
  uint8_t cache;
} task_t;

typedef enum {
  BIC_PERF_MODE_UNRANK = 1 << 0,
  BIC_PERF_MODE_RANK = 1 << 1,
  BIC_PERF_MODE_ALL = BIC_PERF_MODE_UNRANK | BIC_PERF_MODE_RANK,
} bic_perf_mode_t;

static const cli_option_def_t opt_k0 = {
    .opt = {"k0", required_argument, 0, 1000},
    .description = "Starting value for the number of parts.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_k1 = {
    .opt = {"k1", required_argument, 0, 1001},
    .description = "Ending value for the number of parts, inclusive.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_perf_order = {
    .opt = {"order", required_argument, 0, 'o'},
    .description = "Combinatorial orders to test (input is a bitmask).",
    .type_description = "ord",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_perf_unrank_alg = {
    .opt = {"unrank_alg", required_argument, 0, 'a'},
    .description = "Specific unranking strategy for `colex` order to test "
                   "(input is a bitmask).",
    .type_description = "alg",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_perf_rank_alg = {
    .opt = {"rank_alg", required_argument, 0, 'R'},
    .description = "Specific ranking strategy to test (input is a bitmask).",
    .type_description = "alg",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_perf_cache = {
    .opt = {"cache", required_argument, 0, 'c'},
    .description = "Pre-computation of mathematical functions to test (input "
                   "is a bitmask).",
    .type_description = "mode",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_perf_plot = {
    .opt = {"perf", required_argument, 0, 'p'},
    .description = "Measure only ranking, unranking or both functions (input "
                   "is a bitmask).",
    .type_description = "mode",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

size_t generate_tasks(task_t *tasks, const uint16_t orders,
                      const uint16_t unrank_algs, const uint16_t rank_algs,
                      const uint16_t caches, const bic_ctx_t ctx) {
  size_t count = 0;
  for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
    if (!IS_BIT_SET(orders, i) || bic_ctx_set_order(i, ctx)) {
      continue;
    }
    for (uint8_t u = 0; u < BIC_UNRANK_ALG_LENGTH; ++u) {
      if (!IS_BIT_SET(unrank_algs, u) || bic_ctx_set_unrank_alg(u, ctx)) {
        continue;
      }
      for (uint8_t r = 0; r < BIC_RANK_ALG_LENGTH; ++r) {
        if (!IS_BIT_SET(rank_algs, r) || bic_ctx_set_rank_alg(r, ctx)) {
          continue;
        }
        for (uint8_t c = 0; c < BIC_CACHE_LENGTH; ++c) {
          if (!IS_BIT_SET(caches, c) || bic_ctx_set_cache(c, ctx)) {
            continue;
          }
          if (tasks != NULL) {
            task_t task = {
                .order = i, .unrank_alg = u, .rank_alg = r, .cache = c};
            tasks[count] = task;
          }
          ++count;
        }
      }
    }
  }
  return count;
}

void print_perf_header(const task_t *tasks, const size_t count,
                       const uint16_t mode) {
  printf("k");
  for (size_t i = 0; i < count; ++i) {
    const char *order = bic_order_names[tasks[i].order];
    const char *unrank_alg = bic_unrank_alg_names[tasks[i].unrank_alg];
    const char *rank_alg = bic_rank_alg_names[tasks[i].rank_alg];
    const char *cache = bic_cache_names[tasks[i].cache];

    if (mode & BIC_PERF_MODE_UNRANK) {
      printf(" %s-unrank:%s-%s-ucycles", order, unrank_alg, cache);
    }

    if (mode & BIC_PERF_MODE_RANK) {
      printf(" %s-rank:%s-%s-rcycles", order, rank_alg, cache);
    }
  }
  printf("\n");
}

int32_t main(int32_t argc, char **argv) {
  bic_error_t error = BIC_SUCCESS;
  uint32_t m = 0;
  uint32_t k0 = 0;
  uint32_t k1 = 0;
  uint16_t orders = 0;
  uint16_t unrank_algs = 0;
  uint16_t rank_algs = 0;
  uint16_t caches = 0;
  uint16_t mode = 0;
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_target, &m},
      {&opt_k0, &k0},
      {&opt_k1, &k1},
      {&opt_perf_order, &orders},
      {&opt_perf_unrank_alg, &unrank_algs},
      {&opt_perf_rank_alg, &rank_algs},
      {&opt_perf_cache, &caches},
      {&opt_perf_plot, &mode},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  error = parse_args(argc, argv, options);
  if (error != BIC_SUCCESS) {
    return error;
  }

  if ((orders & ~((1u << BIC_ORDER_LENGTH) - 1)) ||
      (unrank_algs & ~((1u << BIC_UNRANK_ALG_LENGTH) - 1)) ||
      (caches & ~((1u << BIC_CACHE_LENGTH) - 1)) ||
      (mode & ~BIC_PERF_MODE_ALL) ||
      (rank_algs & ~((1u << BIC_RANK_ALG_LENGTH) - 1))) {
    return BIC_ERROR_ARGUMENT;
  }

  srandom(seed);

  if (orders == 0) {
    orders = (1 << BIC_ORDER_LENGTH) - 1;
  }
  if (unrank_algs == 0) {
    unrank_algs = (1 << BIC_UNRANK_ALG_LENGTH) - 1;
  }
  if (rank_algs == 0) {
    rank_algs = (1 << BIC_RANK_ALG_LENGTH) - 1;
  }
  if (caches == 0) {
    caches = (1 << BIC_CACHE_LENGTH) - 1;
  }
  if (mode == 0) {
    mode = BIC_PERF_MODE_ALL;
  }

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return BIC_ERROR_ALLOCATION;
  }
  size_t count =
      generate_tasks(NULL, orders, unrank_algs, rank_algs, caches, ctx);
  task_t *tasks = (task_t *)malloc(count * sizeof(task_t));

  if (tasks == NULL) {
    error = BIC_ERROR_ALLOCATION;
    goto cleanup;
  }

  if (m == 0 || k0 == 0 || k0 > k1 || k1 > INT32_MAX || iterations == 0 ||
      count == 0) {
    print_help(argv[0], options);
    error = BIC_ERROR_ARGUMENT;
    goto cleanup;
  }

  generate_tasks(tasks, orders, unrank_algs, rank_algs, caches, ctx);

  bic_ctx_set_strategy(BIC_STRATEGY_GEN, ctx);

  print_perf_header(tasks, count, mode);

  for (uint32_t k = k0; k <= k1; ++k) {
    printf("%u", k);
    for (size_t i = 0; i < count; ++i) {
      bic_ctx_set_order(tasks[i].order, ctx);
      bic_ctx_set_unrank_alg(tasks[i].unrank_alg, ctx);
      bic_ctx_set_rank_alg(tasks[i].rank_alg, ctx);
      bic_ctx_set_cache(tasks[i].cache, ctx);

      uint32_t n = 0;
      uint32_t d = 0;
      bic_run_strategy(m, &n, k, &d, ctx);
      const bic_error_t prepared = n == 0 || n == UINT32_MAX
                                       ? BIC_ERROR_RANGE
                                       : bic_precompute(n, k, d, ctx);
      if (prepared != BIC_SUCCESS) {
        if (mode & BIC_PERF_MODE_UNRANK) {
          printf(" nan");
        }
        if (mode & BIC_PERF_MODE_RANK) {
          printf(" nan");
        }
        bic_free_precomputed(ctx);
        error = prepared;
        continue;
      }

      long double ucycles = 0;
      long double rcycles = 0;

      const bic_error_t measured = bic_run_round_trips(
          n, k, d, iterations, NULL, &ucycles, NULL, &rcycles, ctx);
      if (measured != BIC_SUCCESS) {
        error = measured;
        goto cleanup;
      }

      if (mode & BIC_PERF_MODE_UNRANK) {
        printf(" %.2Lf", ucycles / iterations);
      }

      if (mode & BIC_PERF_MODE_RANK) {
        printf(" %.2Lf", rcycles / iterations);
      }

      bic_free_precomputed(ctx);
    }
    printf("\n");
  }

cleanup:
  free(tasks);
  bic_ctx_destroy(ctx);

  return error;
}
