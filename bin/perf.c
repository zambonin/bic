#include <stdio.h>

#include "extra.h"

#define IS_BIT_SET(mask, idx) ((mask) & (1 << (idx)))

typedef struct {
  uint8_t order;
  uint8_t algorithm;
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
    .type_description = "uint16_t",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_k1 = {
    .opt = {"k1", required_argument, 0, 1001},
    .description = "Ending value for the number of parts, inclusive.",
    .type_description = "uint16_t",
    .parser = parse_uint16,
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

static const cli_option_def_t opt_perf_alg = {
    .opt = {"algorithm", required_argument, 0, 'a'},
    .description = "Specific unranking strategy for `colex` order to test "
                   "(input is a bitmask).",
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

size_t count_tasks(const uint16_t orders, const uint16_t algorithms,
                   const uint16_t caches) {
  size_t count = 0;
  for (size_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
    if (!IS_BIT_SET(orders, i)) {
      continue;
    }
    const uint16_t a =
        (i == BIC_ORDER_COLEX) ? __builtin_popcount(algorithms) : 1;
    count += a * __builtin_popcount(caches);
  }
  return count;
}

void generate_tasks(task_t *tasks, const uint16_t orders,
                    const uint16_t algorithms, const uint16_t caches) {
  size_t count = 0;
  for (uint8_t i = 0; i < BIC_ORDER_LENGTH; ++i) {
    if (!IS_BIT_SET(orders, i)) {
      continue;
    }

    const uint16_t a =
        (i == BIC_ORDER_COLEX) ? algorithms : (1 << BIC_ALG_DEFAULT);
    for (uint8_t j = 0; j < BIC_ALG_LENGTH; ++j) {
      if (!IS_BIT_SET(a, j)) {
        continue;
      }

      for (uint8_t k = 0; k < BIC_CACHE_LENGTH; ++k) {
        if (!IS_BIT_SET(caches, k)) {
          continue;
        }

        task_t t = {.order = i, .algorithm = j, .cache = k};
        tasks[count++] = t;
      }
    }
  }
}

void print_perf_header(const task_t *tasks, const size_t count,
                       const uint16_t mode) {
  printf("k");
  for (size_t i = 0; i < count; ++i) {
    const char *order = bic_order_names[tasks[i].order];
    const char *algorithm = bic_unrank_alg_names[tasks[i].algorithm];
    const char *cache = bic_cache_names[tasks[i].cache];

    if (mode & BIC_PERF_MODE_UNRANK) {
      printf(" %s-%s-%s-ucycles", order, algorithm, cache);
    }

    if (mode & BIC_PERF_MODE_RANK) {
      printf(" %s-%s-%s-rcycles", order, algorithm, cache);
    }
  }
  printf("\n");
}

int32_t main(int32_t argc, char **argv) {
  bool error = false;
  uint16_t m = 0;
  uint16_t k0 = 0;
  uint16_t k1 = 0;
  uint16_t orders = 0;
  uint16_t algorithms = 0;
  uint16_t caches = 0;
  uint16_t mode = 0;
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
      {&opt_target, &m},
      {&opt_k0, &k0},
      {&opt_k1, &k1},
      {&opt_perf_order, &orders},
      {&opt_perf_alg, &algorithms},
      {&opt_perf_cache, &caches},
      {&opt_perf_plot, &mode},
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);

  srandom(seed);

  if (orders == 0) {
    orders = (1 << BIC_ORDER_LENGTH) - 1;
  }
  if (algorithms == 0) {
    algorithms = (1 << BIC_ALG_LENGTH) - 1;
  }
  if (caches == 0) {
    caches = (1 << BIC_CACHE_LENGTH) - 1;
  }
  if (mode == 0) {
    mode = BIC_PERF_MODE_ALL;
  }

  size_t count = count_tasks(orders, algorithms, caches);

  bic_ctx_t ctx = bic_ctx_init();
  task_t *tasks = (task_t *)malloc(count * sizeof(task_t));

  if (m == 0 || k0 == 0 || k1 == 0 || count == 0 || ctx == NULL ||
      tasks == NULL) {
    print_help(argv[0], options);
    error = true;
    goto cleanup;
  }

  generate_tasks(tasks, orders, algorithms, caches);

  bic_ctx_set_strategy(BIC_STRATEGY_GEN, ctx);

  print_perf_header(tasks, count, mode);

  for (uint16_t k = k0; k <= k1; ++k) {
    printf("%u", k);
    for (size_t i = 0; i < count; ++i) {
      bic_ctx_set_order(tasks[i].order, ctx);
      bic_ctx_set_unrank_alg(tasks[i].algorithm, ctx);
      bic_ctx_set_cache(tasks[i].cache, ctx);

      uint16_t n = 0;
      uint16_t d = 0;
      bic_run_strategy(m, &n, k, &d, ctx);
      bic_precompute(n, k, d, ctx);

      long double ucycles = 0;
      long double rcycles = 0;

      bic_run_round_trips(n, k, d, iterations, NULL, &ucycles, NULL, &rcycles,
                          ctx);

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
