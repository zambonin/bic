#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "api.h"
#include "extra.h"

#define INVALID_PARAM                                                          \
  if (fprintf(stderr, "Invalid parameter.\n")) {                               \
    return 1;                                                                  \
  }

static const struct option long_options[] = {
    {"sum", required_argument, 0, 'n'},
    {"parts", required_argument, 0, 'k'},
    {"bound", required_argument, 0, 'd'},
    {"order", required_argument, 0, 'o'},
    {"algorithm", required_argument, 0, 'a'},
    {"iterations", required_argument, 0, 'i'},
    {"cache", required_argument, 0, 'c'},
    {"target", required_argument, 0, 'm'},
    {"strategy", required_argument, 0, 's'},
    {"randomness", required_argument, 0, 'r'},
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
    "  -i, --iterations=<uint32_t>\n"
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
    "  -s, --strategy=<mode>\n"
    "         Use <mode> to choose `n` and `d` according to `m` and `k`.\n"
    "         Available options are:\n"
    "           * `gen` (minimize `d` at all costs);\n"
    "           * `ver` (minimize `n` at all costs).\n"
    "\n"
    "  -r, --randomness=<uint32_t>\n"
    "         Set the seed for the random number generator.\n";

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, uint16_t *m,
                   uint32_t *seed, bic_ctx_t *ctx) {
  for (;;) {
    int c = getopt_long(argc, argv, "n:k:d:o:a:i:c:m:s:r:", long_options, NULL);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'n':
      *n = strtol(optarg, NULL, 0);
      break;
    case 'k':
      *k = strtol(optarg, NULL, 0);
      break;
    case 'd':
      *d = strtol(optarg, NULL, 0);
      break;
    case 'o':
      if (bic_ctx_set_order_by_name(ctx, optarg))
        INVALID_PARAM;
      break;
    case 'a':
      if (bic_ctx_set_unrank_alg_by_name(ctx, optarg))
        INVALID_PARAM;
      break;
    case 'i':
      *iterations = strtol(optarg, NULL, 0);
      break;
    case 'c':
      if (bic_ctx_set_cache_by_name(ctx, optarg))
        INVALID_PARAM;
      break;
    case 'm':
      *m = strtol(optarg, NULL, 0);
      break;
    case 's':
      if (bic_ctx_set_strategy_by_name(ctx, optarg))
        INVALID_PARAM;
      break;
    case 'r':
      *seed = strtol(optarg, NULL, 0);
      break;
    default:
      return 1;
    }
  }

  if (*m != 0 && *k != 0 && *n == 0 && *d == 0) {
    bic_run_strategy(ctx, *m, n, *k, d);
  } else if (*n == 0 || *k == 0 || *d == 0) {
    if (fprintf(stderr, help_text, argv[0])) {
      return 1;
    }
  } else if (*k * *d < *n) {
    INVALID_PARAM;
  }

  return 0;
}

int32_t main(int32_t argc, char **argv) {
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint16_t m = 0;
  uint32_t iterations = 1;
  uint32_t seed = (uint32_t)time(NULL);

  bic_ctx_t *ctx = bic_ctx_init();

  if (parse_args(argc, argv, &n, &k, &d, &iterations, &m, &seed, ctx) > 0) {
    return 1;
  }

  srandom(seed);

  bic_precompute(ctx, n, k, d);

  long double utime = 0;
  long double ucycles = 0;

  long double rtime = 0;
  long double rcycles = 0;

  bic_run_round_trips(ctx, n, k, d, iterations, &utime, &ucycles, &rtime,
                      &rcycles);

  bic_print_config(ctx, n, k, d);
  bic_print_perf_stats(iterations, utime, ucycles, rtime, rcycles);

  bic_ctx_destroy(ctx);

  return 0;
}
