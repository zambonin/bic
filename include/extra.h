#ifndef EXTRA_H
#define EXTRA_H

#include <getopt.h>
#include <time.h>

#include "api.h"
#include "math.h"

static const uint32_t NS_TO_SEC = 1000000000;

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

uint64_t cycles(void);

uint16_t bits_fit_bic(const uint16_t n, const uint16_t k, const uint16_t d);

uintx random_rank(const uint16_t n, const uint16_t k, const uint16_t d);

uintx import_from_unsigned_char(const uint8_t *bytes, const size_t len);

void bic_run_round_trips(const uint16_t n, const uint16_t k, const uint16_t d,
                         const uint32_t iterations, long double *utime,
                         long double *ucycles, long double *rtime,
                         long double *rcycles, const bic_ctx_t ctx);

void bic_run_single_round_trip(const uint16_t n, const uint16_t k,
                               const uint16_t d, const bic_ctx_t ctx);

void bic_print_config(const uint16_t n, const uint16_t k, const uint16_t d,
                      const bic_ctx_t ctx);

void bic_print_perf_stats(const uint32_t it, const long double utime,
                          const long double ucycles, const long double rtime,
                          const long double rcycles);

typedef struct cli_option_def_s {
  struct option opt;
  const char *description;
  const char *type_description;
  uint8_t (*parser)(void *value, const char *arg);
  const char **sub_option_names;
  const char **sub_option_infos;
  size_t num_sub_options;
} cli_option_def_t;

typedef struct {
  const cli_option_def_t *def;
  void *value;
} cli_option_t;

void parse_args(int32_t argc, char **argv, cli_option_t options[]);

uint8_t parse_uint16(void *value, const char *arg);

uint8_t parse_uint32(void *value, const char *arg);

uint8_t parse_string(void *value, const char *arg);

void print_help(const char *name, cli_option_t options[]);

static const cli_option_def_t opt_sum = {
    .opt = {"sum", required_argument, 0, 'n'},
    .description = "Target sum.",
    .type_description = "uint16_t",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_parts = {
    .opt = {"parts", required_argument, 0, 'k'},
    .description = "Number of parts.",
    .type_description = "uint16_t",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const cli_option_def_t opt_bound = {
    .opt = {"bound", required_argument, 0, 'd'},
    .description = "Upper bound of each part, inclusive.",
    .type_description = "uint16_t",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const char *bic_order_sub_options[BIC_ORDER_LENGTH] = {
    "co-lexicographic order",
    "strong minimal-change Gray order",
    "recursive block order due to Miracle-Yilek",
    "spiral order centered on mean",
};

static const cli_option_def_t opt_order = {
    .opt = {"order", required_argument, 0, 'o'},
    .description = "Combinatorial order of choice.",
    .type_description = "ord",
    .parser = parse_string,
    .sub_option_names = bic_order_names,
    .sub_option_infos = bic_order_sub_options,
    .num_sub_options = BIC_ORDER_LENGTH,
};

static const char *bic_unrank_alg_sub_options[BIC_ALG_LENGTH] = {
    "linear search over #C(n, k, d) calculated on demand",
    "linear search reusing partial sums of a single #C(n, k, d)",
    "linear search over a pre-calculated array of accumulated sums of #C(n, k, "
    "d)",
    "binary search over a pre-calculated array of accumulated sums of #C(n, k, "
    "d)",
    "binary search over accumulated sums of #C(n, k, d) calculated directly on "
    "demand",
};

static const cli_option_def_t opt_algorithm = {
    .opt = {"algorithm", required_argument, 0, 'a'},
    .description = "Specific unranking strategy for `colex` order.",
    .type_description = "alg",
    .parser = parse_string,
    .sub_option_names = bic_unrank_alg_names,
    .sub_option_infos = bic_unrank_alg_sub_options,
    .num_sub_options = BIC_ALG_LENGTH,
};

static const cli_option_def_t opt_iterations = {
    .opt = {"iterations", required_argument, 0, 'i'},
    .description = "Number of iterations.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const char *bic_cache_sub_options[BIC_CACHE_LENGTH] = {
    "no pre-computation enabled",
    "binomial coefficients",
    "#C(n, k, d) for intermediate parameters",
    "#C(n, k, d) for intermediate parameters (gaussian)",
    "accumulated sums of #C(n, k, d)",
    "accumulated sums of #C(n, k, d) (gaussian)",
};

static const cli_option_def_t opt_cache = {
    .opt = {"cache", required_argument, 0, 'c'},
    .description = "Pre-computation of mathematical functions.",
    .type_description = "mode",
    .parser = parse_string,
    .sub_option_names = bic_cache_names,
    .sub_option_infos = bic_cache_sub_options,
    .num_sub_options = BIC_CACHE_LENGTH,
};

static const cli_option_def_t opt_target = {
    .opt = {"target", required_argument, 0, 'm'},
    .description = "Target security level, in bits.",
    .type_description = "uint16_t",
    .parser = parse_uint16,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

static const char *bic_strategy_sub_options[BIC_STRATEGY_LENGTH] = {
    "minimize `d` at all costs",
    "minimize `n` at all costs",
    "small pseudorandom parameters",
};

static const cli_option_def_t opt_strategy = {
    .opt = {"strategy", required_argument, 0, 's'},
    .description = "Method to choose `n` and `d` according to `m` and `k`.",
    .type_description = "strat",
    .parser = parse_string,
    .sub_option_names = bic_strategy_names,
    .sub_option_infos = bic_strategy_sub_options,
    .num_sub_options = BIC_STRATEGY_LENGTH,
};

static const cli_option_def_t opt_randomness = {
    .opt = {"randomness", required_argument, 0, 'r'},
    .description = "Set the seed for the random number generator.",
    .type_description = "uint32_t",
    .parser = parse_uint32,
    .sub_option_names = NULL,
    .sub_option_infos = NULL,
    .num_sub_options = 0,
};

#endif // EXTRA_H
