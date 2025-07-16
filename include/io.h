#ifndef IO_H
#define IO_H

#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "common.h"

#define INVALID_PARAM                                                          \
  if (fprintf(stderr, "Invalid parameter.\n")) {                               \
    return 1;                                                                  \
  }

enum {
  PRINT_STATS = 6,
  PRINT_ACCESS = 7,
  PRINT_BUILD = 8,
};

extern int print_type;
extern uint32_t *access_pattern;
extern uint32_t access_pattern_rows;
extern uint32_t access_pattern_cols;

void pprint(const uint16_t n, const uint16_t k, const uint16_t d,
            const uint32_t it, const long double utime,
            const long double ucycles, const long double rtime,
            const long double rcycles);

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, order *ord, uint16_t *m,
                   strategy_func *strategy);

#endif
