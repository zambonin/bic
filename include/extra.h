#ifndef EXTRA_H
#define EXTRA_H

#include "api.h"

void bic_run_round_trips(bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                         const uint16_t d, const uint32_t iterations,
                         long double *utime, long double *ucycles,
                         long double *rtime, long double *rcycles);

void bic_run_single_round_trip(bic_ctx_t *ctx, const uint16_t n,
                               const uint16_t k, const uint16_t d);

void bic_print_config(const bic_ctx_t *ctx, const uint16_t n, const uint16_t k,
                      const uint16_t d);

void bic_print_perf_stats(const uint32_t it, const long double utime,
                          const long double ucycles, const long double rtime,
                          const long double rcycles);

#endif // EXTRA_H
