#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "extra.h"
#include "utils.h"

#define PERF(total_time, total_cycles, logic, var)                             \
  struct timespec var##_tstart;                                                \
  struct timespec var##_tstop;                                                 \
                                                                               \
  clock_gettime(CLOCK_MONOTONIC_RAW, &var##_tstart);                           \
  uint64_t var##_cstart = cycles();                                            \
                                                                               \
  logic;                                                                       \
                                                                               \
  uint64_t var##_cstop = cycles();                                             \
  clock_gettime(CLOCK_MONOTONIC_RAW, &var##_tstop);                            \
                                                                               \
  total_time +=                                                                \
      (long double)(var##_tstop.tv_sec - var##_tstart.tv_sec) * NS_TO_SEC +    \
      (long double)(var##_tstop.tv_nsec - var##_tstart.tv_nsec);               \
  total_cycles += var##_cstop - var##_cstart;

// https://github.com/sphincs/sphincsplus/blob/7ec789ac/ref/test/cycles.c
uint64_t cycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

void check_valid_bounded_composition(const uint32_t *c, const uint16_t n,
                                     const uint16_t k, const uint16_t d) {
  uint32_t sum = 0;
  for (uint16_t i = 0; i < k; ++i) {
    assert(c[i] <= d);
    sum += c[i];
  }
  assert(sum == n);
}

void bic_run_round_trips(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                         const uint16_t d, const uint32_t iterations,
                         long double *utime, long double *ucycles,
                         long double *rtime, long double *rcycles) {
  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));

  for (uint32_t it = 0; it < iterations; ++it) {
    memset(comp, 0, k * sizeof(uint32_t));

    const uintx r = random_rank(n, k, d);

    if (utime && ucycles) {
      PERF(*utime, *ucycles, bic_unrank(ctx, comp, n, k, d, r), unrank);
    } else {
      bic_unrank(ctx, comp, n, k, d, r);
    }

    uintx rr;
    if (rtime && rcycles) {
      PERF(*rtime, *rcycles, rr = bic_rank(ctx, n, k, d, comp), rank);
    } else {
      rr = bic_rank(ctx, n, k, d, comp);
    }

    assert(r == rr);

    check_valid_bounded_composition(comp, n, k, d);
  }

  free(comp);
}

void bic_run_single_round_trip(bic_ctx_t *ctx, const uint16_t n,
                               const uint16_t k, const uint16_t d) {
  bic_run_round_trips(ctx, n, k, d, 1, NULL, NULL, NULL, NULL);
}

void bic_print_config(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                      const uint16_t d) {
  const char *order_name = bic_order_find_ctx.names[ctx->order];
  const char *unrank_alg_name = bic_unrank_alg_find_ctx.names[ctx->unrank_alg];
  const char *strategy_name = bic_strategy_find_ctx.names[ctx->strategy];
  printf("n=%-5u k=%-5u d=%-5u m=%-5d b=%-5.0f c=%-2u o=%-5s a=%-7s s=%-6s\n",
         n, k, d, bits_fit_bic(n, k, d), BIT_LENGTH, ctx->cache_type,
         order_name, unrank_alg_name, strategy_name);
}

void bic_print_perf_stats(const uint32_t it, const long double utime,
                          const long double ucycles, const long double rtime,
                          const long double rcycles) {
  printf("i=%-5u, unrank avg = %14.2Lf ns, %14.2Lf cyc., rank avg = %14.2Lf "
         "ns, %14.2Lf cyc.\n",
         it, utime / it, ucycles / it, rtime / it, rcycles / it);
}
