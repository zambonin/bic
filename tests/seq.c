#include <stdio.h>
#include <string.h>

#include "extra.h"

static const uint32_t LIMIT = 10;

int32_t main(int32_t argc, char **argv) {
  uint32_t iterations = 8;
  uint32_t seed = (uint32_t)time(NULL);

  cli_option_t options[] = {
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
          for (uint32_t it = 0; it < iterations; ++it) {
            for (uint32_t n = 0; n < LIMIT; ++n) {
              for (uint32_t k = 0; k < LIMIT; ++k) {
                for (uint32_t d = 0; d < LIMIT; ++d) {
                  if (k * d < n) {
                    continue;
                  }
                  if (bic_precompute(n, k, d, ctx) != 0) {
                    continue;
                  }
                  bic_run_round_trips(n, k, d, 1, NULL, NULL, NULL, NULL, ctx);
                  bic_free_precomputed(ctx);
                }
              }
            }
          }
        }
      }
    }
  }

  bic_ctx_destroy(ctx);

  puts("sequential round trips passed");
  return 0;
}
