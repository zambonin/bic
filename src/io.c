#include "io.h"

int print_type = PRINT_STATS;
uint32_t *access_pattern;
uint32_t access_pattern_rows = 0;
uint32_t access_pattern_cols = 0;

void pprint(const uint16_t n, const uint16_t k, const uint16_t d,
            const uint32_t iterations, const long double total_time,
            const long double total_cycles) {
  if (print_type == PRINT_STATS) {
    printf("n = %5d, k = %5d, d = %5d, i = %5u, m = %10.4Lf, "
           "avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",
           n, k, d, iterations, lg_bic(n, k, d), total_time / iterations,
           total_cycles / iterations);
  } else if (print_type == PRINT_ACCESS) {
    PRINT_MATRIX_INT(access_pattern, access_pattern_rows, access_pattern_cols);
  }
}

void free_cache() {
  if (cache_type >= BIN_CACHE) {
    free(bin_cache);
  }
  if (cache_type >= COMB_CACHE) {
    free(comb_cache);
  }
  if (cache_type >= ACC_COMB_CACHE) {
    for (uint16_t i = 0; i < acc_cache_rows; ++i) {
      for (uint16_t j = 0; j < acc_cache_cols; ++j) {
        free(GET_CACHE(acc_cache, i, j));
      }
    }
    free(acc_cache);
  }
  free(access_pattern);
}

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, unrank_func *unrank,
                   rank_func *rank, uint16_t *m, strategy_func *strategy) {
  for (;;) {
    int c = getopt_long(argc, argv, "n:k:d:a:i:c:m:p:s:", long_options, NULL);
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
    case 'a':
      if (strcmp(optarg, "colex") == 0) {
        *unrank = colex;
        *rank = colex_rank;
      } else if (strcmp(optarg, "colexbs") == 0) {
        *unrank = colex_bs;
        *rank = colex_rank;
      } else if (strcmp(optarg, "colexpart") == 0) {
        *unrank = colex_part;
        *rank = colex_rank;
      } else if (strcmp(optarg, "colexdbcs") == 0) {
        *unrank = colex_dbcs;
        *rank = colex_rank;
      } else if (strcmp(optarg, "gray") == 0) {
        *unrank = gray;
        *rank = gray_rank;
      } else if (strcmp(optarg, "rbo") == 0) {
        *unrank = rbo;
        *rank = rbo_rank;
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
