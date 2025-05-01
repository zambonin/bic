#ifndef IO_H
#define IO_H

#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "common.h"

#define PRINT_MATRIX_GENERIC(expr, n_rows, n_cols)                             \
  for (uint32_t row = 0; row < n_rows; ++row) {                                \
    for (uint32_t col = 0; col < n_cols; ++col) {                              \
      expr;                                                                    \
    }                                                                          \
    printf("\n");                                                              \
  }

#define PRINT_MATRIX_LG(var, n_rows, n_cols)                                   \
  PRINT_MATRIX_GENERIC(printf("%8.2Lf ", lg(GET_CACHE(var, row, col))),        \
                       n_rows, n_cols)

#define PRINT_MATRIX_INT(var, n_rows, n_cols)                                  \
  PRINT_MATRIX_GENERIC(printf("%8d ", GET_CACHE(var, row, col)), n_rows, n_cols)

#define PRINT_CACHE_STATS(var, type, desc)                                     \
  if (cache_type == type && print_type == PRINT_ACCESS) {                      \
    access_pattern = (uint32_t *)calloc(n_rows * n_cols, sizeof(uint32_t));    \
    assert(access_pattern != NULL);                                            \
    access_pattern_rows = n_rows;                                              \
    access_pattern_cols = n_cols;                                              \
  } else if (print_type == PRINT_BUILD) {                                      \
    printf("n = %5d, k = %5d, d = %5d, c = %6s, "                              \
           "avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",                    \
           n, k, d, desc, total_time, total_cycles);                           \
  } else if (cache_type == type && print_type == PRINT_LENGTH) {               \
    PRINT_MATRIX_LG(var, n_rows, n_cols);                                      \
  }

#define INVALID_PARAM                                                          \
  if (fprintf(stderr, "Invalid parameter.\n")) {                               \
    return 1;                                                                  \
  }

enum {
  PRINT_STATS = 4,
  PRINT_ACCESS = 5,
  PRINT_LENGTH = 6,
  PRINT_BUILD = 7,
};

static const struct option long_options[] = {
    {"sum", required_argument, 0, 'n'},
    {"parts", required_argument, 0, 'k'},
    {"bound", required_argument, 0, 'd'},
    {"order", required_argument, 0, 'o'},
    {"algorithm", required_argument, 0, 'a'},
    {"iterations", required_argument, 0, 'i'},
    {"cache", required_argument, 0, 'c'},
    {"target", required_argument, 0, 'm'},
    {"print", required_argument, 0, 'p'},
    {"strategy", required_argument, 0, 's'},
    {0, 0, 0, 0},
};

static const char *help_text =
    "Ranking and unranking algorithms for integer compositions with summands\n"
    "bounded in number and size, referred to as C(n, k, d).\n"
    "\n"
    "Usage: %s [OPTIONS]\n"
    "  -n, --sum=<uint16_t>\n"
    "         Target sum.\n"
    "\n"
    "  -k, --parts=<uint16_t>\n"
    "         Number of parts.\n"
    "\n"
    "  -d, --bound=<uint16_t>\n"
    "         Upper bound of each part, inclusive.\n"
    "\n"
    "  -o, --order=<ord>\n"
    "         Use <ord> as the combinatorial order of choice.\n"
    "         Available options are:\n"
    "           * `colex` (co-lexicographic order);\n"
    "           * `gray` (strong minimal-change Gray order);\n"
    "           * `rbo` (recursive block order due to Miracle-Yilek).\n"
    "\n"
    "  -a, --algorithm=<alg>\n"
    "         Use <alg> as the specific unranking strategy for `colex`.\n"
    "         Available options are:\n"
    "           * `default` (linear search over #C(n, k, d) calculated on\n"
    "               demand);\n"
    "           * `ps` (linear search reusing partial sums of a single\n"
    "               #C(n, k, d));\n"
    "           * `al` (linear search over a pre-calculated array of\n"
    "               accumulated sums of #C(n, k, d));\n"
    "           * `ab` (binary search over a pre-calculated array of\n"
    "               accumulated sums of #C(n, k, d));\n"
    "           * `ad` (binary search over accumulated sums of #C(n, k, d)\n"
    "               calculated directly on demand).\n"
    "\n"
    "  -i, --iterations=<uint64_t>\n"
    "         Number to repeatedly unrank random integers.\n"
    "\n"
    "  -c, --cache=<mode>\n"
    "         Use <mode> to pre-compute mathematical functions.\n"
    "         Available options are:\n"
    "           * `bin` (binomial coefficients);\n"
    "           * `comb` (#C(n, k, d) for intermediate parameters);\n"
    "           * `acc` (accumulated sums of #C(n, k, d)).\n"
    "\n"
    "  -m, --target=<uint16_t>\n"
    "         Target security level, in bits.\n"
    "\n"
    "  -p, --print=<mode>\n"
    "         Show information about the pre-computed cache.\n"
    "         Available options are:\n"
    "           * `access` (count how frequently elements are accessed);\n"
    "           * `length` (show bit length of all elements);\n"
    "           * `build` (measure performance of building the cache).\n"
    "\n"
    "  -s, --strategy=<mode>\n"
    "         Use <mode> to choose `n` and `d` according to `m` and `k`.\n"
    "         Available options are:\n"
    "           * `gen` (minimize `d` at all costs);\n"
    "           * `ver` (minimize `n` at all costs).\n";

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
