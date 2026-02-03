#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "extra.h"
#include "strat.h"

static const uint32_t SMALL = 64;
static const uint32_t LIMIT = 1000;

uint32_t uniform_1_to_64(void) { return (random() % SMALL) + 1; }

void gen_small_params_mingen(uint32_t *n, uint32_t *k, uint32_t *d) {
  do {
    uint32_t m = uniform_1_to_64();
    *k = uniform_1_to_64();
    mingen(m, n, *k, d);
  } while (*n > LIMIT || !*n);
}

void gen_small_params_minver(uint32_t *n, uint32_t *k, uint32_t *d) {
  do {
    uint32_t m = uniform_1_to_64();
    *k = uniform_1_to_64();
    minver(m, n, *k, d);
  } while (*n > LIMIT || !*n);
}

void gen_small_params_random(uint32_t *n, uint32_t *k, uint32_t *d) {
  *k = uniform_1_to_64();
  *d = uniform_1_to_64();
  *n = (random() % ((*k * *d + 1) / 2)) + 1;
}

void bic_run_random_params(uint32_t *n, uint32_t *k, uint32_t *d,
                           const bic_ctx_t ctx) {
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

  cli_option_t options[] = {
      {&opt_iterations, &iterations},
      {&opt_randomness, &seed},
      {NULL, NULL},
  };

  parse_args(argc, argv, options);

  srandom(seed);

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }

  for (uint8_t i = 0; i < BIC_ORDER_LENGTH; i++) {
    bic_ctx_set_order(i, ctx);
    for (uint8_t j = 0; j < BIC_UNRANK_ALG_LENGTH; j++) {
      if (bic_ctx_set_unrank_alg(j, ctx)) {
        continue;
      }
      for (uint8_t r = 0; r < BIC_RANK_ALG_LENGTH; r++) {
        if (bic_ctx_set_rank_alg(r, ctx)) {
          continue;
        }
        for (uint8_t c = 0; c < BIC_CACHE_LENGTH; c++) {
          bic_ctx_set_cache(c, ctx);
          for (uint8_t s = 0; s < BIC_STRATEGY_LENGTH; s++) {
            bic_ctx_set_strategy(s, ctx);
            for (uint32_t it = 0; it < iterations; ++it) {
              uint32_t n = 0;
              uint32_t k = 0;
              uint32_t d = 0;
              bic_run_random_params(&n, &k, &d, ctx);
              bic_print_config(n, k, d, ctx);
              printf("\n");
              bic_precompute(n, k, d, ctx);
              bic_run_single_round_trip(n, k, d, ctx);
              bic_free_precomputed(ctx);
            }
          }
        }
      }
    }
  }

  bic_ctx_destroy(ctx);
  return 0;
}
