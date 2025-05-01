#include "io.h"
#include "cache.h"
#include "colex.h"
#include "gray.h"
#include "math.h"
#include "rbo.h"

int print_type = PRINT_STATS;
uint32_t *access_pattern;
uint32_t access_pattern_rows = 0;
uint32_t access_pattern_cols = 0;

void pprint(const uint16_t n, const uint16_t k, const uint16_t d,
            const uint32_t it, const long double utime,
            const long double ucycles, const long double rtime,
            const long double rcycles) {
  if (print_type == PRINT_STATS) {
    printf("n = %5d, k = %5d, d = %5d, i = %5u, m = %10.4Lf, "
           "unrank avg = %14.2Lf ns, %14.2Lf cyc., "
           "rank avg = %14.2Lf ns, %14.2Lf cyc.\n",
           n, k, d, it, lg_bic(n, k, d), utime / it, ucycles / it, rtime / it,
           rcycles / it);
  } else if (print_type == PRINT_ACCESS) {
    PRINT_MATRIX_INT(access_pattern, access_pattern_rows, access_pattern_cols);
  }
}

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, order *ord, uint16_t *m,
                   strategy_func *strategy) {
  for (;;) {
    int c = getopt_long(argc, argv, "n:k:d:o:a:i:c:m:p:s:", long_options, NULL);
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
      } else if (strcmp(optarg, "comb") == 0) {
        cache_type = COMB_CACHE;
      } else if (strcmp(optarg, "acc") == 0) {
        cache_type = ACC_COMB_CACHE;
      } else
        INVALID_PARAM;
      break;
    case 'm':
      *m = strtol(optarg, NULL, 0);
      break;
    case 'p':
      if (strcmp(optarg, "none") == 0) {
        print_type = PRINT_STATS;
      } else if (strcmp(optarg, "access") == 0) {
        print_type = PRINT_ACCESS;
      } else if (strcmp(optarg, "length") == 0) {
        print_type = PRINT_LENGTH;
      } else if (strcmp(optarg, "build") == 0) {
        print_type = PRINT_BUILD;
      } else
        INVALID_PARAM;
      break;
    case 's':
      if (strcmp(optarg, "gen") == 0) {
        *strategy = mingen;
      } else if (strcmp(optarg, "ver") == 0) {
        *strategy = minver;
      } else
        INVALID_PARAM;
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
