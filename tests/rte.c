#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "extra.h"

typedef struct {
  const char *name;
  const uint32_t n;
  const uint32_t k;
  const uint32_t d;
  const bic_cache_t cache;
} param_ctx_t;

static const uint8_t KEY_LENGTH = 32;

static const param_ctx_t USE_CASES[] = {
    {"10x10", 32, 100, 255, BIC_CACHE_SMALL_COMB},
    {"16x16", 24, 256, 255, BIC_CACHE_NONE},
    {"32x32", 16, 1024, 255, BIC_CACHE_NONE},
    {"Exams", 24, 300, 100, BIC_CACHE_NONE},
    {"Salaries", 64, 30, 100000, BIC_CACHE_SMALL_COMB},
    {"Ratings", 8, 5000, 4, BIC_CACHE_NONE},
    {NULL, 0, 0, 0, BIC_CACHE_NONE},
};

static void xxor(uintx *rop, const uintx msg, const uint8_t *const key,
                 const size_t len) {
  *rop = (msg ^ import_from_unsigned_char(key, len));
}

static int32_t
run_rte_round_trips(const bic_ctx_t ctx, const uint32_t n, const uint32_t k,
                    const uint32_t d, const uint32_t iterations,
                    long double *const etime, long double *const ecycles,
                    long double *const dtime, long double *const dcycles) {
  uint32_t *comp = (uint32_t *)calloc(k, sizeof(uint32_t));
  uint32_t *ret = (uint32_t *)calloc(k, sizeof(uint32_t));
  if (comp == NULL || ret == NULL) {
    free(ret);
    free(comp);
    return 1;
  }

  const uint32_t bits_len = bits_fit_bic(n, k, d);
  const uint32_t len = (bits_len + 7) / 8;
  uint8_t key[KEY_LENGTH];
  uint8_t mkey[KEY_LENGTH];
  for (size_t i = 0; i < KEY_LENGTH; ++i) {
    key[i] = (uint8_t)random();
  }

  const uint32_t rem = bits_len % 8;
  const uint32_t mask = (rem == 0) ? 0xFF : (1U << rem) - 1U;
  for (size_t i = 0; i < len; ++i) {
    mkey[i] = key[i];
  }
  mkey[0] &= mask;

  int32_t status = 0;
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

    if (r != rank || rank != dec) {
      status = 1;
      goto cleanup;
    }
    for (size_t i = 0; i < k; ++i) {
      if (comp[i] != ret[i]) {
        status = 1;
        goto cleanup;
      }
      comp[i] = ret[i] = 0;
    }
  }

cleanup:
  free(ret);
  free(comp);
  return status;
}

static int32_t run_test_combination(bic_order_t order, bic_cache_t cache,
                                    uint32_t n, uint32_t k, uint32_t d,
                                    uint32_t iterations) {
  if (!bic_params_supported(n, k, d, cache) ||
      bits_fit_bic(n, k, d) > KEY_LENGTH * 8) {
    fprintf(stderr, "unsupported RTE parameters: n=%u, k=%u, d=%u\n", n, k, d);
    return 1;
  }

  bic_ctx_t ctx = bic_ctx_init();
  if (ctx == NULL) {
    return 1;
  }
  int32_t status = 1;
  long double etime = 0, dtime = 0;
  long double ecycles = 0, dcycles = 0;
  if (bic_ctx_set_order(order, ctx) != 0 ||
      bic_ctx_set_cache(cache, ctx) != 0) {
    goto cleanup;
  }
  if (bic_precompute(n, k, d, ctx) != 0 &&
      (errno != ERANGE || bic_ctx_set_cache(BIC_CACHE_NONE, ctx) != 0 ||
       bic_precompute(n, k, d, ctx) != 0)) {
    fprintf(stderr, "RTE precomputation failed: n=%u, k=%u, d=%u\n", n, k, d);
    goto cleanup;
  }
  if (ctx->cache_type != cache) {
    fprintf(stderr,
            "RTE: unsupported statistical cache; using no cache for "
            "n=%u, k=%u, d=%u\n",
            n, k, d);
  }

  bic_print_config(n, k, d, ctx);

  if (run_rte_round_trips(ctx, n, k, d, iterations, &etime, &ecycles, &dtime,
                          &dcycles) != 0) {
    fputs("RTE round trip failed\n", stderr);
    goto cleanup;
  }

  bic_print_perf_stats(iterations, etime, ecycles, dtime, dcycles);
  status = 0;

cleanup:
  bic_free_precomputed(ctx);
  bic_ctx_destroy(ctx);
  return status;
}

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

  if (iterations == 0) {
    fputs("RTE iterations must be positive\n", stderr);
    return 1;
  }

  for (uint8_t i = 0; USE_CASES[i].name != NULL; ++i) {
    const uint32_t n = USE_CASES[i].n;
    const uint32_t k = USE_CASES[i].k;
    const uint32_t d = USE_CASES[i].d;

    printf("Case: %s (n=%u, k=%u, d=%u)\n", USE_CASES[i].name, n, k, d);
    if (run_test_combination(BIC_ORDER_COLEX, USE_CASES[i].cache, n, k, d,
                             iterations) != 0) {
      return 1;
    }
    printf("\n");
    if (run_test_combination(BIC_ORDER_SPIRAL, USE_CASES[i].cache, n, k, d,
                             iterations) != 0) {
      return 1;
    }
    printf("\n");
  }

  return 0;
}
