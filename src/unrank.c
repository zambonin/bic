// NOLINTBEGIN(readability-identifier-length, llvmlibc-restrict-system-libc-headers, concurrency-mt-unsafe, hicpp-no-assembler, misc-include-cleaner, clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling, cppcoreguidelines-avoid-non-const-global-variables, bugprone-easily-swappable-parameters)

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum {
  NO_CACHE = 0,
  BIN_CACHE = 1,
  COMB_CACHE = 2,
  ACC_COMB_CACHE = 3,
  BIT_LENGTH = 512,
  NS_TO_SEC = 1000000000,
};

typedef unsigned _BitInt(BIT_LENGTH) uintx;
typedef _BitInt(BIT_LENGTH) intx;
typedef void (*unrank_func)(uint16_t *, const uint16_t, const uint16_t,
                            const uint16_t, uintx);

void lex_it(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void lex_plus_it(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void gray_it(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void enup_it(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);

static int cache_type = NO_CACHE;
static uint16_t param_k = 0;
static uintx *comb_cache;
static uintx *bin_cache;

static const struct option long_options[] = {
    {"sum", required_argument, 0, 'n'},
    {"parts", required_argument, 0, 'k'},
    {"bound", required_argument, 0, 'd'},
    {"algorithm", required_argument, 0, 'a'},
    {"iterations", required_argument, 0, 'i'},
    {"cache", required_argument, 0, 'c'},
    {0, 0, 0, 0},
};

static uint64_t cpucycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

uint32_t min(const uint32_t a, const uint32_t b) { return (a < b) ? a : b; }

uintx mfac_uiui(uint64_t n) {
  uintx x = n + (n == 0);
  while (n > 2) {
    x *= --n;
  }
  return x;
}

uintx bin_uiui(uint64_t n, uint64_t k) {
  uintx rr = k <= n;

  if (k > (n >> 1U)) {
    k = (k <= n) ? n - k : 0;
  }

  uintx kfac = mfac_uiui(k);

  for (; k > 0; --k) {
    rr *= n--;
  }

  return rr / kfac;
}

uintx bin(const uint64_t n, const uint64_t k) {
  if (cache_type == BIN_CACHE) {
    return bin_cache[n * param_k + k];
  }

  return bin_uiui(n, k);
}

uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx inner = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = bin(k, i);
    right = bin(n - (d + 1) * i + k - 1, k - 1);
    inner = left * right;
    if (i & 1U) {
      inner = -inner;
    }
    rop += inner;

    if (partial_sums != NULL) {
      partial_sums[i] = inner;
    }
  }

  return (uintx)rop;
}

uintx inner_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return inner_bic_with_sums(n, k, d, NULL);
}

uintx bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  if (cache_type == COMB_CACHE) {
    (void)d;
    return comb_cache[k + param_k * n];
  }

  return inner_bic(n, k, d);
}

void build_bin_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  param_k = k;
  (void)d;
  bin_cache = calloc((uint32_t)((n + k + 1) * param_k), sizeof(uintx));
  assert(bin_cache != NULL);

  for (uint32_t row = 0; row < (uint32_t) (n + k + 1); ++row) {
    bin_cache[0 + k * row] = 1;
    bin_cache[1 + k * row] = row;

    for (uint16_t col = 2; col < k; ++col) {
      bin_cache[col + k * row] = bin_uiui(row, col);
    }
  }
}

void build_comb_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  param_k = k + 1;
  comb_cache = calloc((uint32_t)((n + 1) * param_k), sizeof(uintx));
  assert(comb_cache != NULL);

  comb_cache[0] = 1;
  for (uint16_t row = 1; row <= k; ++row) {
    for (uint16_t col = 0; col <= n; ++col) {
      comb_cache[row + param_k * col] = inner_bic(col, row, d);
    }
  }
}

void lex_it(uint16_t *rop, const uint16_t n, const uint16_t k, const uint16_t d,
            uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0;
         count = bic(it_n - part, i, d), part < min(it_n, d) && rank >= count;
         ++part, rank -= count) {
    }
  }

  rop[0] = it_n;
}

void lex_plus_it(uint16_t *rop, const uint16_t n, const uint16_t k,
                 const uint16_t d, uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;

  uint16_t j = min(k - 1, it_n / (d + 1));
  intx *prev_sum = calloc(j + 1, sizeof(intx));

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    uintx left = 0;
    intx right = inner_bic_with_sums(it_n, i, d, prev_sum);

    for (part = 0; part < min(it_n, d) && rank >= (uintx) right; ++part) {
      left = right;
      uint16_t numerator = (it_n - part);
      uint16_t denominator = (it_n - part) + i - 1;

      j = min(i, (it_n - part - 1) / (d + 1));
      for (uint16_t m = 0; m <= j; ++m) {
        prev_sum[m] *= numerator;
        prev_sum[m] /= denominator;
        right += prev_sum[m];
        numerator -= (d + 1);
        denominator -= (d + 1);
      }
    }

    rank -= left;
  }

  rop[0] = it_n;

  free(prev_sum);
}

void gray_it(uint16_t *rop, const uint16_t n, const uint16_t k,
             const uint16_t d, uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0;
         count = bic(it_n - part, i, d), part < min(it_n, d) && rank >= count;
         ++part, rank -= count) {
    }
    if (part & 1U) {
      rank = count - 1 - rank;
    }
  }

  rop[0] = it_n;
}

void enup_it(uint16_t *rop, const uint16_t n, const uint16_t k,
             const uint16_t d, uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;
  uintx count = 0;

  rank = bic(it_n, k, d) - 1 - rank;
  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0;
         count = bic(it_n - part, i, d), part < min(it_n, d) && rank >= count;
         ++part, rank -= count) {
    }
    if (!(part & 1U)) {
      rank = count - 1 - rank;
    }
  }

  rop[0] = it_n;
}

void check_valid_bounded_composition(const uint16_t *c, const uint16_t n,
                                     const uint16_t k, const uint16_t d) {
  uint16_t sum = 0;
  for (uint16_t i = 0; i < k; ++i) {
    assert(c[i] <= d);
    sum += c[i];
  }
  assert(sum == n);
}

int32_t parse_args(int32_t argc, char **argv, uint16_t *n, uint16_t *k,
                   uint16_t *d, uint32_t *iterations, unrank_func *unrank) {
  for (;;) {
    int c = getopt_long(argc, argv, "n:k:d:a:i:c:", long_options, NULL);
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
      if (strcmp(optarg, "lex") == 0) {
        *unrank = lex_it;
      } else if (strcmp(optarg, "lexplus") == 0) {
        *unrank = lex_plus_it;
      } else if (strcmp(optarg, "gray") == 0) {
        *unrank = gray_it;
      } else if (strcmp(optarg, "enup") == 0) {
        *unrank = enup_it;
      } else if (fprintf(stderr, "Invalid parameter.\n")) {
        return 1;
      }
      break;
    case 'i':
      *iterations = strtol(optarg, NULL, 0);
      break;
    case 'c':
      if (strcmp(optarg, "bin") == 0) {
        cache_type = BIN_CACHE;
      } else if (strcmp(optarg, "comb") == 0) {
        cache_type = COMB_CACHE;
      } else if (strcmp(optarg, "acc") == 0) {
        cache_type = ACC_COMB_CACHE;
      } else if (fprintf(stderr, "Invalid parameter.\n")) {
        return 1;
      }
      break;
    default:
      return 1;
    }
  }

  if (*n == 0 || *k == 0 || *d == 0) {
    int32_t ret = fprintf(
        stderr,
        "Unranking algorithms for integer compositions with summands\n"
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
        "  -a, --algorithm=<alg>\n"
        "         Use <alg> as the unranking algorithm of choice.\n"
        "         Available options are:\n"
        "           * `lex` (lexicographic order);\n"
        "           * `lexplus` (lexicographic order reusing partial sums);\n"
        "           * `gray` (minimal change order, or Gray code);\n"
        "           * `enup` (even numbers up, odd numbers down).\n"
        "\n"
        "  -i, --iterations=<uint64_t>\n"
        "         Number to repeatedly unrank random integers.\n"
        "\n"
        "  -c, --cache=<mode>\n"
        "         Use <mode> to pre-compute mathematical functions.\n"
        "         Available options are:\n"
        "           * `bin` (binomial coefficients);\n"
        "           * `comb` (#C(n, k, d) for intermediate parameters);\n"
        "           * `acc` (accumulated sums of #C(n, k, d)).\n",
        argv[0]);
    (void)ret;
    return 1;
  }

  if (*k * *d < *n) {
    if (fprintf(stderr, "Invalid parameter.\n")) {
      return 1;
    }
  }

  return 0;
}

int32_t main(int32_t argc, char **argv) {
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint32_t iterations = 1;
  unrank_func unrank = lex_it;

  if (parse_args(argc, argv, &n, &k, &d, &iterations, &unrank) > 0) {
    return 1;
  }

  if (cache_type == BIN_CACHE) {
    build_bin_cache(n, k, d);
  } else if (cache_type == COMB_CACHE) {
    build_comb_cache(n, k, d);
  }

  long double comb_lg = logl(inner_bic(n, k, d)) / log(2);
  uint16_t len = (1 + comb_lg) / sizeof(uint64_t);

  long double total_time = 0;
  struct timespec time_start;
  struct timespec time_stop;

  long double total_cycles = 0;
  uint64_t cycle_start = 0;
  uint64_t cycle_stop = 0;

  srandom(time(NULL));
  for (uint32_t it = 0; it < iterations; ++it) {
    uint8_t *message = calloc(len, sizeof(uint8_t));
    for (uint16_t i = 0; i < len; ++i) {
      message[i] = random();
    }

    uintx rank = 0;
    unsigned char *ptr = (unsigned char *)&rank;

    for (uint16_t i = 0; i < len; ++i) {
      ptr[len - 1 - i] = message[i];
    }
    uint16_t *comp = calloc(k, sizeof(uint16_t));

    clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);
    cycle_start = cpucycles();

    (*unrank)(comp, n, k, d, rank);

    cycle_stop = cpucycles();
    clock_gettime(CLOCK_MONOTONIC_RAW, &time_stop);

    total_time +=
        (long double)(time_stop.tv_sec - time_start.tv_sec) * NS_TO_SEC +
        (long double)(time_stop.tv_nsec - time_start.tv_nsec);
    total_cycles += cycle_stop - cycle_start;

    check_valid_bounded_composition(comp, n, k, d);

    free(comp);
    free(message);
  }

  printf("n = %5d, k = %5d, d = %5d, i = %5u, m = %10.4Lf, "
         "avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",
         n, k, d, iterations, comb_lg, total_time / iterations,
         total_cycles / iterations);

  free(comb_cache);
  free(bin_cache);

  return 0;
}

// NOLINTEND(readability-identifier-length, llvmlibc-restrict-system-libc-headers, concurrency-mt-unsafe, hicpp-no-assembler, misc-include-cleaner, clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling, cppcoreguidelines-avoid-non-const-global-variables, bugprone-easily-swappable-parameters)
