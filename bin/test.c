#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "api.h"
#include "extra.h"
#include "strat.h"

static const uint16_t SMALL = 64;
static const uint16_t LIMIT = 1000;

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

int32_t parse_args(int32_t argc, char **argv, uint32_t *iterations,
                   uint32_t *seed) {
  for (;;) {
    int c = getopt_long(argc, argv, "i:s:", long_options, NULL);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      *iterations = strtol(optarg, NULL, 10);
      break;
    case 's':
      *seed = strtol(optarg, NULL, 10);
      break;
    default:
      fprintf(stderr, help_text, argv[0]);
      return 1;
    }
  }

  return 0;
}

uint16_t uniform_1_to_64(void) { return (random() % SMALL) + 1; }

void gen_small_params_mingen(uint16_t *n, uint16_t *k, uint16_t *d) {
  do {
    uint16_t m = uniform_1_to_64();
    *k = uniform_1_to_64();
    mingen(m, n, *k, d);
  } while (*n > LIMIT || !*n);
}

void gen_small_params_minver(uint16_t *n, uint16_t *k, uint16_t *d) {
  do {
    uint16_t m = uniform_1_to_64();
    *k = uniform_1_to_64();
    minver(m, n, *k, d);
  } while (*n > LIMIT || !*n);
}

void gen_small_params_random(uint16_t *n, uint16_t *k, uint16_t *d) {
  *k = uniform_1_to_64();
  *d = uniform_1_to_64();
  *n = (random() % ((*k * *d + 1) / 2)) + 1;
}

void bic_run_random_params(bic_ctx_t *ctx, uint16_t *n, uint16_t *k,
                           uint16_t *d) {
  switch (ctx->strategy) {
  case BIC_STRATEGY_GEN:
    gen_small_params_mingen(n, k, d);
    break;
  case BIC_STRATEGY_VER:
    gen_small_params_minver(n, k, d);
    break;
  case BIC_STRATEGY_RANDOM:
    gen_small_params_random(n, k, d);
    break;
  default:
    return;
  }
}

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  if (parse_args(argc, argv, &iterations, &seed) > 0) {
    return 1;
  }

  srandom(seed);

  bic_ctx_t *ctx = bic_ctx_init();

  for (uint8_t i = 0; i < BIC_ORDER_LENGTH; i++) {
    bic_ctx_set_order(ctx, i);
    for (uint8_t j = 0; j < BIC_ALG_LENGTH; j++) {
      if (bic_ctx_set_unrank_alg(ctx, j)) {
        continue;
      }
      for (uint8_t c = 0; c < BIC_CACHE_LENGTH; c++) {
        bic_ctx_set_cache(ctx, c);
        for (uint8_t s = 0; s < BIC_STRATEGY_LENGTH; s++) {
          bic_ctx_set_strategy(ctx, s);
          for (uint32_t it = 0; it < iterations; ++it) {
            uint16_t n = 0, d = 0, k = 0;
            bic_run_random_params(ctx, &n, &k, &d);
            bic_print_config(ctx, n, k, d);
            bic_precompute(ctx, n, k, d);
            bic_run_single_round_trip(ctx, n, k, d);
            bic_free_precomputed(ctx);
          }
        }
      }
    }
  }

  bic_ctx_destroy(ctx);

  return 0;
}
