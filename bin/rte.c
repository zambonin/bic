#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "extra.h"

typedef struct {
  const char *name;
  const uint32_t k;
  const uint32_t d;
} param_ctx_t;

static const uint8_t KEY_LENGTH = 32;

static const param_ctx_t USE_CASES[] = {
    {"10x10", 100, 255},
    // {"16x16", 256, 255},
    // {"32x32", 1024, 255},
    // {"Exams", 300, 100},
    // {"Salaries", 30, 100000},
    // {"Ratings", 5000, 4},
    {NULL, 0, 0},
};

void xxor(uintx *rop, const uintx msg, const uint8_t *const key,
          const size_t len) {
  *rop = (msg ^ import_from_unsigned_char(key, len));
}

void run_rte_round_trips(const bic_ctx_t ctx, const uint32_t n,
                         const uint32_t k, const uint32_t d,
                         const uint32_t iterations, long double *const etime,
                         long double *const ecycles, long double *const dtime,
                         long double *const dcycles) {
  uint32_t *comp = (uint32_t *)calloc(k, sizeof(uint32_t));
  uint32_t *ret = (uint32_t *)calloc(k, sizeof(uint32_t));
  uint32_t bits_len = bits_fit_bic(n, k, d);
  uint32_t len = (bits_len + 7) / 8;

  uint8_t *key = (uint8_t *)calloc(KEY_LENGTH, sizeof(uint8_t));
  for (size_t i = 0; i < KEY_LENGTH; ++i) {
    key[i] = random();
  }

  uint8_t *mkey = (uint8_t *)calloc(len, sizeof(uint8_t));
  uint32_t rem = bits_len % 8;
  uint32_t mask = (rem == 0) ? 0xFF : (1 << rem) - 1;
  for (size_t i = 0; i < len; ++i) {
    mkey[i] = key[i % KEY_LENGTH];
  }
  mkey[0] &= mask;

  for (uint32_t it = 0; it < iterations; ++it) {
    const uintx r = random_rank(n, k, d);
    bic_unrank(comp, n, k, d, r, ctx);

    uintx enc = 0;
    uintx rank = 0;
    PERF(*etime, *ecycles, rank = bic_rank(n, k, d, comp, ctx);
         xxor(&enc, rank, mkey, len), encipher);

    uintx dec = 0;
    PERF(*dtime, *dcycles, xxor(&dec, enc, mkey, len);
         bic_unrank(ret, n, k, d, dec, ctx), decipher);

    assert(r == rank && rank == dec);
    for (size_t i = 0; i < k; ++i) {
      assert(comp[i] == ret[i]);
      comp[i] = ret[i] = 0;
    }
  }

  free(mkey);
  free(key);
  free(ret);
  free(comp);
}

void run_test_combination(bic_order_t order, bic_cache_t cache, uint32_t n,
                          uint32_t k, uint32_t d, uint32_t iterations) {
  bic_ctx_t ctx = bic_ctx_init();
  bic_ctx_set_order(order, ctx);
  bic_ctx_set_cache(cache, ctx);

  bic_print_config(n, k, d, ctx);
  bic_precompute(n, k, d, ctx);

  long double etime = 0, dtime = 0;
  long double ecycles = 0, dcycles = 0;
  run_rte_round_trips(ctx, n, k, d, iterations, &etime, &ecycles, &dtime,
                      &dcycles);

  bic_print_perf_stats(iterations, etime, ecycles, dtime, dcycles);

  bic_free_precomputed(ctx);
  bic_ctx_destroy(ctx);
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

  for (uint8_t i = 0; USE_CASES[i].name != NULL; ++i) {
    const uint32_t k = USE_CASES[i].k;
    const uint32_t d = USE_CASES[i].d;
    const uint32_t np = ((uint32_t)k * d + 1) / 2;
    const uint32_t n = (np > UINT32_MAX) ? UINT32_MAX : np;

    run_test_combination(BIC_ORDER_RBO, BIC_CACHE_SMALL_COMB, n, k, d,
                         iterations);
  }

  bic_ctx_destroy(ctx);

  return 0;
}