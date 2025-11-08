#include <stdio.h>
#include <string.h>

#include "extra.h"
#include "strat.h"

static const uint16_t SMALL = 64;
static const uint16_t LIMIT = 1000;

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

void bic_run_random_params(uint16_t *n, uint16_t *k, uint16_t *d,
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
    for (uint8_t j = 0; j < BIC_ALG_LENGTH; j++) {
      if (bic_ctx_set_unrank_alg(j, ctx)) {
        continue;
      }
      for (uint8_t c = 0; c < BIC_CACHE_LENGTH; c++) {
        bic_ctx_set_cache(c, ctx);
        for (uint8_t s = 0; s < BIC_STRATEGY_LENGTH; s++) {
          bic_ctx_set_strategy(s, ctx);
          for (uint32_t it = 0; it < iterations; ++it) {
            uint16_t n = 0;
            uint16_t k = 0;
            uint16_t d = 0;
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

  bic_ctx_destroy(ctx);

  return 0;
}
