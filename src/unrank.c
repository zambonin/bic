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
  PRINT_STATS = 4,
  PRINT_ACCESS = 5,
  PRINT_LENGTH = 6,
  PRINT_BUILD = 7,
  BIT_LENGTH = 320,
  NS_TO_SEC = 1000000000,
};

#if defined(BOOST_FIX_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::checked_uint512_t;
using intx = boost::multiprecision::checked_int512_t;
#elif defined(BOOST_ARB_INT)
#include <boost/multiprecision/cpp_int.hpp>

using uintx = boost::multiprecision::cpp_int;
using intx = boost::multiprecision::cpp_int;
#elif defined(BOOST_MPZ_INT)
#include <boost/multiprecision/gmp.hpp>

using uintx = boost::multiprecision::mpz_int;
using intx = boost::multiprecision::mpz_int;
#elif defined(BOOST_TOM_INT)
#include <boost/multiprecision/tommath.hpp>

using uintx = boost::multiprecision::tom_int;
using intx = boost::multiprecision::tom_int;
#endif

#if defined(BITINT)
typedef unsigned _BitInt(BIT_LENGTH) uintx;
typedef _BitInt(BIT_LENGTH) intx;
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/random.hpp>
#include <random>

static std::random_device rd;
#endif

typedef void (*unrank_func)(uint16_t *, const uint16_t, const uint16_t,
                            const uint16_t, uintx);
typedef uintx (*math_func)(const uint16_t, const uint16_t, const uint16_t);

void colex(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void colex_bs(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void colex_part(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void gray(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);
void rbo(uint16_t *rop, uint16_t n, uint16_t k, uint16_t d, uintx rank);

static int cache_type = NO_CACHE;
static uintx *bin_cache;
static uint16_t bin_cache_cols = 0;
static uintx *comb_cache;
static uint16_t comb_cache_cols = 0;
static uintx **acc_cache;
static uint16_t acc_cache_rows = 0;
static uint16_t acc_cache_cols = 0;

static int print_type = PRINT_STATS;
static uint32_t *access_pattern;
static uint32_t access_pattern_rows = 0;
static uint32_t access_pattern_cols = 0;

// from https://stackoverflow.com/a/40245287
typedef struct {
  const void *key;
  int (*const compar)(const void *, const void *);
  void *last_visited;
} bsearch_insertion_state;

static int bsearch_insertion_compare(const void *a, const void *b) {
  bsearch_insertion_state *state = (bsearch_insertion_state *)a;
  state->last_visited = (void *)b;
  return state->compar(state->key, b);
}

size_t bsearch_insertion(const void *key, const void *base, size_t nel,
                         size_t width,
                         int (*compar)(const void *, const void *)) {
  bsearch_insertion_state state = {key, compar, NULL};
  bsearch(&state, base, nel, width, bsearch_insertion_compare);

  uintx *last = (uintx *)state.last_visited;
  if (*((uintx *)key) < *last) {
    --last;
  }

  return last - (uintx *)base;
}

static int cmp(const void *c, const void *d) {
  uintx a = *(uintx *)c;
  uintx b = *(uintx *)d;

  if (a < b) {
    return -1;
  }

  if (a >= b) {
    return 1;
  }

  return 0;
}

#define GET_CACHE(var, i, j) *(var + ((i) * var##_cols) + (j))

#define GET_CACHE_OR_CALC(type, var, math)                                     \
  if (cache_type >= type) {                                                    \
    if (cache_type == type && print_type == PRINT_ACCESS) {                    \
      (GET_CACHE(access_pattern, n, k))++;                                     \
    }                                                                          \
    return GET_CACHE(var, n, k);                                               \
  }                                                                            \
  return math(n, k, d);

#define PRINT_MATRIX_GENERIC(expr, n_rows, n_cols)                             \
  for (uint32_t row = 0; row < n_rows; ++row) {                                \
    for (uint32_t col = 0; col < n_cols; ++col) {                              \
      expr;                                                                    \
    }                                                                          \
    printf("\n");                                                              \
  }

#define PRINT_MATRIX_LG(var, n_rows, n_cols)                                   \
  PRINT_MATRIX_GENERIC(printf("%8.2Lf ", logl2(GET_CACHE(var, row, col))),     \
                       n_rows, n_cols)

#define PRINT_MATRIX_INT(var, n_rows, n_cols)                                  \
  PRINT_MATRIX_GENERIC(printf("%8d ", GET_CACHE(var, row, col)), n_rows, n_cols)

#define MEASURE(total_time, total_cycles, logic)                               \
  struct timespec time_start;                                                  \
  struct timespec time_stop;                                                   \
                                                                               \
  clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);                             \
  uint64_t cycle_start = cpucycles();                                          \
                                                                               \
  logic;                                                                       \
                                                                               \
  uint64_t cycle_stop = cpucycles();                                           \
  clock_gettime(CLOCK_MONOTONIC_RAW, &time_stop);                              \
                                                                               \
  total_time +=                                                                \
      (long double)(time_stop.tv_sec - time_start.tv_sec) * NS_TO_SEC +        \
      (long double)(time_stop.tv_nsec - time_start.tv_nsec);                   \
  total_cycles += cycle_stop - cycle_start;

#define BUILD_CACHE(rows, cols, var, type, desc, logic)                        \
  uint32_t n_rows = rows;                                                      \
  uint32_t n_cols = cols;                                                      \
  var##_cols = n_cols;                                                         \
                                                                               \
  var = (uintx *)calloc(n_rows * n_cols, sizeof(uintx));                       \
  assert(var != NULL);                                                         \
                                                                               \
  if (cache_type == type && print_type == PRINT_ACCESS) {                      \
    access_pattern = (uint32_t *)calloc(n_rows * n_cols, sizeof(uint32_t));    \
    assert(access_pattern != NULL);                                            \
    access_pattern_rows = n_rows;                                              \
    access_pattern_cols = n_cols;                                              \
  }                                                                            \
                                                                               \
  long double total_time = 0;                                                  \
  long double total_cycles = 0;                                                \
                                                                               \
  MEASURE(total_time, total_cycles, logic);                                    \
                                                                               \
  if (print_type == PRINT_BUILD) {                                             \
    printf("n = %5d, k = %5d, d = %5d, c = %6s, "                              \
           "avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",                    \
           n, k, d, desc, total_time, total_cycles);                           \
  } else if (cache_type == type && print_type == PRINT_LENGTH) {               \
    PRINT_MATRIX_LG(var, n_rows, n_cols);                                      \
  }

static const struct option long_options[] = {
    {"sum", required_argument, 0, 'n'},
    {"parts", required_argument, 0, 'k'},
    {"bound", required_argument, 0, 'd'},
    {"algorithm", required_argument, 0, 'a'},
    {"iterations", required_argument, 0, 'i'},
    {"cache", required_argument, 0, 'c'},
    {"target", required_argument, 0, 'm'},
    {"print", required_argument, 0, 'p'},
    {0, 0, 0, 0},
};

static uint64_t cpucycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

uint32_t min(const uint32_t a, const uint32_t b) { return (a < b) ? a : b; }

#if defined(BITINT)
long double logl2(const uintx u) { return logl(u) / log(2); }
#else
long double logl2(const uintx u) {
  return (long double)log2(boost::multiprecision::cpp_bin_float_100(u));
}
#endif

// from FXT: aux0/binomial.h
uintx bin_uiui(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  if (k > n)
    return 0;
  if ((k == 0) || (k == n))
    return 1;

  uint16_t kk = k;
  if (2 * k > n)
    kk = n - k;

  uintx b = n - kk + 1;
  uintx f = b;
  for (uintx j = 2; j <= kk; ++j) {
    ++f;
    b *= f;
    b /= j;
  }
  return b;
}

uintx bin(const uint16_t n, const uint16_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(BIN_CACHE, bin_cache, bin_uiui);
}

uintx inner_bic_with_sums(const uint16_t n, const uint16_t k, const uint16_t d,
                          intx *partial_sums, math_func bin_impl) {
  if (n == 0) {
    return 1;
  }

  intx rop = 0;
  intx inner = 0;
  uintx left = 0;
  uintx right = 0;

  uint16_t j = min(k, n / (d + 1));
  for (uint16_t i = 0; i <= j; ++i) {
    left = bin_impl(k, i, d);
    right = bin_impl(n - (d + 1) * i + k - 1, k - 1, d);
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
  return inner_bic_with_sums(n, k, d, NULL, bin);
}

uintx bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(COMB_CACHE, comb_cache, inner_bic);
}

uintx *inner_acc(const uint16_t n, const uint16_t k, const uint16_t d) {
  size_t length = d + 3;
  size_t i = 0;
  uintx *rop = (uintx *)calloc(length, sizeof(uintx));
  uintx sum = 0;

  rop[0] = 0;
  for (; i <= min(n, d); ++i) {
    sum += bic(n - i, k, d);
    rop[i + 1] = sum;
  }
  rop[length - 1] = i;

  return rop;
}

uintx *acc(const uint16_t n, const uint16_t k, const uint16_t d) {
  GET_CACHE_OR_CALC(ACC_COMB_CACHE, acc_cache, inner_acc);
}

void build_bin_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  (void)d;
  BUILD_CACHE(n + k + 1, k, bin_cache, BIN_CACHE, "bin", {
    GET_CACHE(bin_cache, 0, 0) = 1;
    for (uint32_t row = 1; row < n_rows; ++row) {
      GET_CACHE(bin_cache, row, 0) = 1;
      for (uint16_t col = 1; col <= min(row, n_cols - 1); ++col) {
        GET_CACHE(bin_cache, row, col) =
            (uintx)GET_CACHE(bin_cache, row - 1, col - 1) +
            (uintx)GET_CACHE(bin_cache, row - 1, col);
      }
    }
  });
}

void build_comb_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  BUILD_CACHE(n + 1, k + 1, comb_cache, COMB_CACHE, "comb", {
    GET_CACHE(comb_cache, 0, 0) = 1;
    for (uint16_t row = 0; row < n_rows; ++row) {
      for (uint16_t col = 1; col < n_cols; ++col) {
        GET_CACHE(comb_cache, row, col) = inner_bic(row, col, d);
      }
    }
  });
}

void build_acc_cache(const uint16_t n, const uint16_t k, const uint16_t d) {
  uint32_t n_rows = n + 1;
  uint32_t n_cols = k;

  acc_cache_rows = n_rows;
  acc_cache_cols = n_cols;

  acc_cache = (uintx **)calloc(n_rows * n_cols, sizeof(uintx *));
  assert(acc_cache != NULL);

  for (uint16_t row = 0; row < n_rows; ++row) {
    for (uint16_t col = 0; col < n_cols; ++col) {
      GET_CACHE(acc_cache, row, col) = inner_acc(row, col, d);
    }
  }
}

long double lg_bic(const uint16_t n, const uint16_t k, const uint16_t d) {
  return logl2(inner_bic_with_sums(n, k, d, NULL, bin_uiui));
}

void get_min_params(const uint16_t m, uint16_t *n, const uint16_t k,
                    uint16_t *d) {
  uint16_t it_d = 0;
  uint16_t it_n = (k * it_d) / 2;

  while (lg_bic(it_n, k, it_d) < m) {
    ++it_d;
    it_n = (k * it_d) / 2;
  }

  while (lg_bic(it_n, k, it_d) >= m) {
    --it_d;
  }
  ++it_d;

  while (lg_bic(it_n, k, it_d) >= m) {
    --it_n;
  }
  ++it_n;

  *n = it_n;
  *d = it_d;
}

void colex(uint16_t *rop, const uint16_t n, const uint16_t k, const uint16_t d,
           uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;
  uintx count = 0;

  if (cache_type == ACC_COMB_CACHE) {
    for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
      uintx *sums = acc(it_n, i, d);
      for (part = 0; count = sums[part + 1], rank >= count; ++part) {
      }
      rank -= sums[part];
    }
  } else {
    for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
      for (part = 0; count = bic(it_n - part, i, d), rank >= count;
           ++part, rank -= count) {
      }
    }
  }

  rop[0] = it_n;
}

void colex_bs(uint16_t *rop, const uint16_t n, const uint16_t k,
              const uint16_t d, uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    uintx *sums = acc(it_n, i, d);
    size_t length = (size_t)sums[d + 2];
    part = bsearch_insertion(&rank, sums, length, sizeof(uintx), cmp);
    rank -= sums[part];
  }

  rop[0] = it_n;
}

void colex_part(uint16_t *rop, const uint16_t n, const uint16_t k,
                const uint16_t d, uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;

  uint16_t j = min(k - 1, it_n / (d + 1));
  intx *prev_sum = (intx *)calloc(j + 1, sizeof(intx));

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    intx left = 0;
    intx right = inner_bic_with_sums(it_n, i, d, prev_sum, bin);

    for (part = 0; rank >= (uintx)right; ++part) {
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

    rank -= (uintx)left;
  }

  rop[0] = it_n;

  free(prev_sum);
}

void gray(uint16_t *rop, const uint16_t n, const uint16_t k, const uint16_t d,
          uintx rank) {
  uint16_t it_n = n;
  uint16_t part = 0;
  uintx count = 0;

  for (uint16_t i = k - 1; i > 0; rop[i] = part, --i, it_n -= part) {
    for (part = 0; count = bic(it_n - part, i, d), rank >= count;
         ++part, rank -= count) {
    }
    if (part & 1U) {
      rank = count - 1 - rank;
    }
  }

  rop[0] = it_n;
}

void inner_rbo(uint16_t *rop, const uint16_t n, const uint16_t k,
               const uint16_t d, uintx rank, uint16_t start) {
  if (k == 1) {
    rop[start] = n;
    return;
  }

  uint16_t left = (uint16_t)(k / 2);
  uint16_t right = k - left;

  uint16_t leftSum = 0;
  uintx partSum = 0;

  uint16_t s = 0;
  for (uintx count = 0; s <= min(n, left * d) && count <= rank;
       ++s, count += bic(n - s + 1, left, d) * bic(s - 1, right, d)) {
    leftSum = s;
    partSum = count;
  }

  uint16_t rightSum = n - leftSum;
  uintx rightPoints = bic(rightSum, right, d);
  uintx leftRank = (rank - partSum) / rightPoints;
  uintx rightRank = (rank - partSum) % rightPoints;

  inner_rbo(rop, leftSum, left, d, leftRank, start);
  inner_rbo(rop, rightSum, right, d, rightRank, start + left);
}

void rbo(uint16_t *rop, const uint16_t n, const uint16_t k, const uint16_t d,
         uintx rank) {
  inner_rbo(rop, n, k, d, rank, 0);
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
                   uint16_t *d, uint32_t *iterations, unrank_func *unrank,
                   uint16_t *m) {
  for (;;) {
    int c = getopt_long(argc, argv, "n:k:d:a:i:c:m:p:", long_options, NULL);
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
      } else if (strcmp(optarg, "colexbs") == 0) {
        *unrank = colex_bs;
      } else if (strcmp(optarg, "colexpart") == 0) {
        *unrank = colex_part;
      } else if (strcmp(optarg, "gray") == 0) {
        *unrank = gray;
      } else if (strcmp(optarg, "rbo") == 0) {
        *unrank = rbo;
      } else if (fprintf(stderr, "Invalid parameter.\n")) {
        return 1;
      }
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
      } else if (fprintf(stderr, "Invalid parameter.\n")) {
        return 1;
      }
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
      } else if (fprintf(stderr, "Invalid parameter.\n")) {
        return 1;
      }
      break;
    default:
      return 1;
    }
  }

  if (*m != 0 && *k != 0 && *n == 0 && *d == 0) {
    get_min_params(*m, n, *k, d);
  } else if (*n == 0 || *k == 0 || *d == 0) {
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
        "           * `colex` (co-lexicographic order);\n"
        "           * `colexpart` (co-lexicographic order reusing partial "
        "sums);\n"
        "           * `gray` (strong minimal-change Gray order).\n"
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
        "           * `build` (measure performance of building the cache).\n",
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

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d) {
#if defined(BITINT)
  long double comb_lg = lg_bic(n, k, d);
  uint16_t len = (1 + comb_lg) / sizeof(uint64_t);
  len += (len == 0);

  uint8_t *message = (uint8_t *)calloc(len, sizeof(uint8_t));
  for (uint16_t i = 0; i < len; ++i) {
    message[i] = random();
  }

  uintx rank = 0;
  unsigned char *ptr = (unsigned char *)&rank;

  for (uint16_t i = 0; i < len; ++i) {
    ptr[len - 1 - i] = message[i];
  }
  free(message);

  return rank % inner_bic(n, k, d);
#else
  boost::random::uniform_int_distribution<uintx> ui(0, inner_bic(n, k, d));
  return ui(rd);
#endif
}

int32_t main(int32_t argc, char **argv) {
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  uint16_t m = 0;
  uint32_t iterations = 1;
  unrank_func unrank = colex;

  if (parse_args(argc, argv, &n, &k, &d, &iterations, &unrank, &m) > 0) {
    return 1;
  }

  srandom(time(NULL));

  if (cache_type >= BIN_CACHE) {
    build_bin_cache(n, k, d);
    if (cache_type >= COMB_CACHE) {
      build_comb_cache(n, k, d);
      if (cache_type >= ACC_COMB_CACHE) {
        build_acc_cache(n, k, d);
      }
    }
  }

  long double total_time = 0;
  long double total_cycles = 0;
  uint16_t *comp = (uint16_t *)malloc(k * sizeof(uint16_t));

  for (uint32_t it = 0; it < iterations; ++it) {
    memset(comp, 0, k * sizeof(uint16_t));

    uintx rank = random_rank(n, k, d);

    MEASURE(total_time, total_cycles, (*unrank)(comp, n, k, d, rank));

    check_valid_bounded_composition(comp, n, k, d);
  }

  if (print_type == PRINT_STATS) {
    printf("n = %5d, k = %5d, d = %5d, i = %5u, m = %10.4Lf, "
           "avg time = %14.2Lf ns, avg cycles = %14.2Lf\n",
           n, k, d, iterations, lg_bic(n, k, d), total_time / iterations,
           total_cycles / iterations);
  } else if (print_type == PRINT_ACCESS) {
    PRINT_MATRIX_INT(access_pattern, access_pattern_rows, access_pattern_cols);
    free(access_pattern);
  }

  if (cache_type >= BIN_CACHE) {
    free(bin_cache);
    if (cache_type >= COMB_CACHE) {
      free(comb_cache);
      if (cache_type >= ACC_COMB_CACHE) {
        for (uint16_t i = 0; i < acc_cache_rows; ++i) {
          for (uint16_t j = 0; j < acc_cache_cols; ++j) {
            free(GET_CACHE(acc_cache, i, j));
          }
        }
        free(acc_cache);
      }
    }
  }

  free(comp);

  return 0;
}
