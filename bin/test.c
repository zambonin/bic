#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cache.h"
#include "colex.h"
#include "common.h"
#include "gray.h"
#include "math.h"
#include "rbo.h"
#include "utils.h"

typedef void (*param_gen_func)(uint16_t *, uint16_t *, uint16_t *);

void gen_params_mingen(uint16_t *n, uint16_t *k, uint16_t *d);
void gen_params_minver(uint16_t *n, uint16_t *k, uint16_t *d);
void gen_params_random(uint16_t *n, uint16_t *k, uint16_t *d);

typedef struct {
  const char *name;
  void (*unrank_func)(uint32_t *, const uint16_t, const uint16_t,
                      const uint16_t, const uintx);
} algo_t;

typedef struct {
  const char *name;
  order ord;
  const algo_t *algos;
  size_t num_algos;
} order_cfg_t;

typedef struct {
  const char *name;
  int strategy;
} cache_cfg_t;

typedef struct {
  const char *name;
  param_gen_func func;
} strategy_cfg_t;

static const algo_t ALGOS_COLEX[] = {{"default", colex_unrank},
                                     {"ps", colex_unrank_part_sums},
                                     {"al", colex_unrank_acc_linear},
                                     {"ab", colex_unrank_acc_bisect},
                                     {"ad", colex_unrank_acc_direct}};

static const order_cfg_t ORDERS[] = {
    {"colex", colex, ALGOS_COLEX, sizeof(ALGOS_COLEX) / sizeof(algo_t)},
    {"gray", gray, NULL, 0},
    {"rbo", rbo, NULL, 0}};

static const cache_cfg_t CACHES[] = {{"none", NO_CACHE},
                                     {"bin", BIN_CACHE},
                                     {"comb", COMB_CACHE},
                                     {"scomb", SMALL_COMB_CACHE},
                                     {"acc", ACC_COMB_CACHE}};

static const strategy_cfg_t STRATEGIES[] = {{"mingen", gen_params_mingen},
                                            {"minver", gen_params_minver},
                                            {"random", gen_params_random}};

static const struct option long_options[] = {
    {"iterations", required_argument, 0, 'i'},
    {"seed", required_argument, 0, 's'},
    {0, 0, 0, 0}};

static const char *help_text =
    "Usage: %s [OPTIONS]\n"
    "  -i, --iterations=<uint32_t>\n"
    "         Set the number of random tests per configuration.\n"
    "\n"
    "  -s, --seed=<uint32_t>\n"
    "         Set the seed for the random number generator.\n";

void gen_params_mingen(uint16_t *n, uint16_t *k, uint16_t *d) {
  do {
    uint16_t m = (random() % 64) + 1;
    *k = (random() % 64) + 1;
    mingen(m, n, *k, d);
  } while (*n > 1000 || !*n);
}

void gen_params_minver(uint16_t *n, uint16_t *k, uint16_t *d) {
  do {
    uint16_t m = (random() % 64) + 1;
    *k = (random() % 64) + 1;
    minver(m, n, *k, d);
  } while (*n > 1000 || !*n);
}

void gen_params_random(uint16_t *n, uint16_t *k, uint16_t *d) {
  *k = (random() % 64) + 1;
  *d = (random() % 64) + 1;
  *n = (random() % ((*k * *d + 1) / 2)) + 1;
}

void run_round_trip(order ord, const uint16_t n, const uint16_t k,
                    const uint16_t d) {
  uintx all = inner_bic_with_sums(n, k, d, NULL, inner_bin);
  if (all == 0) {
    return;
  }

  const uintx r = random_rank(n, k, d);
  uintx rr;

  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));
  assert(comp != NULL);

  (*ord.unrank)(comp, n, k, d, r);
  check_valid_bounded_composition(comp, n, k, d);
  rr = (*ord.rank)(n, k, d, comp);

  assert(r == rr);

  free(comp);
}

void report_test(const char *order_name, const char *algo_name,
                 const char *strat_name, uint32_t n, uint32_t k, uint32_t d) {
  printf("n=%-5u k=%-5u d=%-5u m=%-5d b=%-5.0f c=%-2u o=%-5s a=%-7s s=%-6s\n",
         n, k, d, bits_fit_bic(n, k, d), BIT_LENGTH, cache_type, order_name,
         algo_name, strat_name);
}

void run_suite(uint32_t iterations) {
  const size_t num_orders = sizeof(ORDERS) / sizeof(order_cfg_t);
  const size_t num_caches = sizeof(CACHES) / sizeof(cache_cfg_t);
  const size_t num_strats = sizeof(STRATEGIES) / sizeof(strategy_cfg_t);

  for (size_t i = 0; i < num_orders; ++i) {
    order_cfg_t order_cfg = ORDERS[i];
    size_t num_algos = (order_cfg.algos) ? order_cfg.num_algos : 1;

    for (size_t j = 0; j < num_algos; ++j) {
      order test_order = order_cfg.ord;
      const char *algo_name = "default";

      if (order_cfg.algos) {
        test_order.unrank = order_cfg.algos[j].unrank_func;
        algo_name = order_cfg.algos[j].name;
      }

      for (size_t c = 0; c < num_caches; ++c) {
        cache_type = CACHES[c].strategy;

        for (size_t s = 0; s < num_strats; ++s) {
          strategy_cfg_t strat_cfg = STRATEGIES[s];

          for (uint32_t t = 0; t < iterations; ++t) {
            uint16_t n = 0, d = 0, k = 0;
            strat_cfg.func(&n, &k, &d);
            report_test(order_cfg.name, algo_name, strat_cfg.name, n, k, d);
            build_caches(n, k, d);
            run_round_trip(test_order, n, k, d);
            free_caches();
          }
        }
      }
    }
  }
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  for (;;) {
    int c = getopt_long(argc, argv, "i:s:", long_options, NULL);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      iterations = strtoul(optarg, NULL, 10);
      break;
    case 's':
      seed = strtoul(optarg, NULL, 10);
      break;
    default:
      fprintf(stderr, help_text, argv[0]);
      return 1;
    }
  }

  srandom(seed);
  run_suite(iterations);

  return 0;
}
