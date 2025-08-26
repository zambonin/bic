#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "cache.h"
#include "colex.h"
#include "gray.h"
#include "math.h"
#include "rbo.h"
#include "utils.h"

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

void pprint(const uint16_t n, const uint16_t k, const uint16_t d,
            const uint32_t it, const long double utime,
            const long double ucycles, const long double rtime,
            const long double rcycles) {
  printf("n = %5d, k = %5d, d = %5d, i = %5u, m = %10.4Lf, b = %5f, c = %2u, "
         "unrank avg = %14.2Lf ns, %14.2Lf cyc., "
         "rank avg = %14.2Lf ns, %14.2Lf cyc.\n",
         n, k, d, it, bits_fit_bic(n, k, d), BIT_LENGTH, cache_type, utime / it,
         ucycles / it, rtime / it, rcycles / it);
}

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, order *ord, uint16_t *m,
                   strategy_func *strategy, uint32_t *seed) {
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
      if (strcmp(optarg, "colex") == 0) {
        *ord = colex;
      } else if (strcmp(optarg, "gray") == 0) {
        *ord = gray;
      } else if (strcmp(optarg, "rbo") == 0) {
        *ord = rbo;
      } else
        INVALID_PARAM;
      break;
    case 'a':
      if (strcmp(optarg, "default") == 0) {
        (void)optarg;
      } else if ((*ord).rank != colex_rank) {
        INVALID_PARAM;
      } else if (strcmp(optarg, "al") == 0) {
        (*ord).unrank = colex_unrank_acc_linear;
      } else if (strcmp(optarg, "ab") == 0) {
        (*ord).unrank = colex_unrank_acc_bisect;
      } else if (strcmp(optarg, "ps") == 0) {
        (*ord).unrank = colex_unrank_part_sums;
      } else if (strcmp(optarg, "ad") == 0) {
        (*ord).unrank = colex_unrank_acc_direct;
      } else
        INVALID_PARAM;
      break;
    case 'i':
      *iterations = strtol(optarg, NULL, 0);
      break;
    case 'c':
      if (strcmp(optarg, "none") == 0) {
        cache_type = NO_CACHE;
      } else if (strcmp(optarg, "bin") == 0) {
        cache_type = BIN_CACHE;
      } else if (strcmp(optarg, "rbo") == 0) {
        cache_type = RBO_COMB_CACHE;
      } else if (strcmp(optarg, "comb") == 0) {
        cache_type = COMB_CACHE;
      } else if (strcmp(optarg, "acc") == 0) {
        cache_type = ACC_COMB_CACHE;
      } else if (strcmp(optarg, "scomb") == 0) {
        cache_type = SMALL_COMB_CACHE;
      } else
        INVALID_PARAM;
      break;
    case 'm':
      *m = strtol(optarg, NULL, 0);
      break;
    case 's':
      if (strcmp(optarg, "gen") == 0) {
        *strategy = mingen;
      } else if (strcmp(optarg, "ver") == 0) {
        *strategy = minver;
      } else
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
    (*strategy)(*m, n, *k, d);
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
  uint32_t seed = time(NULL);
  order ord = colex;
  strategy_func strategy = mingen;

  if (parse_args(argc, argv, &n, &k, &d, &iterations, &ord, &m, &strategy,
                 &seed) > 0) {
    return 1;
  }

  srandom(seed);

  build_caches(n, k, d);

  long double utime = 0;
  long double ucycles = 0;

  long double rtime = 0;
  long double rcycles = 0;

  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));

  for (uint32_t it = 0; it < iterations; ++it) {
    memset(comp, 0, k * sizeof(uint32_t));

    const uintx r = random_rank(n, k, d);

    PERF(utime, ucycles, (*ord.unrank)(comp, n, k, d, r), unrank);

    PERF(rtime, rcycles, const uintx rr = (*ord.rank)(n, k, d, comp), rank);

    assert(r == rr);

    check_valid_bounded_composition(comp, n, k, d);
  }

  pprint(n, k, d, iterations, utime, ucycles, rtime, rcycles);

  free_caches();

  free(comp);

  return 0;
}
