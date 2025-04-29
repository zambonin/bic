#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <time.h>

#include "common.h"

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

// from https://stackoverflow.com/a/40245287
typedef struct {
  const void *key;
  void *last_visited;
} bsearch_insertion_state;

uint64_t cycles(void);

uint32_t min(const uint32_t a, const uint32_t b);

int32_t max(const int32_t a, const int32_t b);

void check_valid_bounded_composition(const uint32_t *c, const uint16_t n,
                                     const uint16_t k, const uint16_t d);

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width);

#endif
