#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "extra.h"

uint64_t cycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

void check_valid_bounded_composition(const uint32_t *c, const uint32_t n,
                                     const uint32_t k, const uint32_t d) {
  uint32_t sum = 0;
  for (uint32_t i = 0; i < k; ++i) {
    assert(c[i] <= d);
    sum += c[i];
  }
  assert(sum == n);
}

uint16_t bits_fit_bic(const uint32_t n, const uint32_t k, const uint32_t d) {
  return 1 + ((uint16_t)lg(compute_bic_no_cache(n, k, d)));
}

uintx random_rank(const uint32_t n, const uint32_t k, const uint32_t d) {
  uint16_t len = bits_fit_bic(n, k, d) / sizeof(uint64_t);
  len += (len == 0);

  uint8_t *message = (uint8_t *)calloc(len, sizeof(uint8_t));
  for (uint16_t i = 0; i < len; ++i) {
    message[i] = random();
  }

  uintx rank = import_from_unsigned_char(message, len);

  free(message);
  return rank % compute_bic_no_cache(n, k, d);
}

uintx import_from_unsigned_char(const uint8_t *bytes, const size_t len) {
  uintx rop = 0;

#if defined(BITINT)
  unsigned char *ptr = (unsigned char *)&rop;
  for (uint16_t i = 0; i < len; ++i) {
    ptr[len - 1 - i] = bytes[i];
  }
#elif defined(BOOST_FIX_INT) || defined(BOOST_ARB_INT)
  boost::multiprecision::detail::import_bits_fast(rop, bytes, bytes + len);
#elif defined(BOOST_MPZ_INT)
  mpz_import(rop.backend().data(), len, 1, sizeof(uint8_t), 0, 0, bytes);
#elif defined(BOOST_TOM_INT)
  mp_err err =
      mp_unpack(&rop.backend().data(), len, 1, sizeof(uint8_t), 0, 0, bytes);
  (void)err;
#endif

  return rop;
}

void bic_run_round_trips(const uint32_t n, const uint32_t k, const uint32_t d,
                         const uint32_t iterations, long double *utime,
                         long double *ucycles, long double *rtime,
                         long double *rcycles, const bic_ctx_t ctx) {
  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));
  if (comp == NULL) {
    return;
  }

  for (uint32_t it = 0; it < iterations; ++it) {
    memset(comp, 0, k * sizeof(uint32_t));

    const uintx r = random_rank(n, k, d);

    if (utime || ucycles) {
      long double dummy_time;
      long double dummy_cycles;
      long double *t = utime ? utime : &dummy_time;
      long double *c = ucycles ? ucycles : &dummy_cycles;
      PERF(*t, *c, bic_unrank(comp, n, k, d, r, ctx), unrank);
    } else {
      bic_unrank(comp, n, k, d, r, ctx);
    }

    uintx rr;
    if (rtime || rcycles) {
      long double dummy_time;
      long double dummy_cycles;
      long double *t = rtime ? rtime : &dummy_time;
      long double *c = rcycles ? rcycles : &dummy_cycles;
      PERF(*t, *c, rr = bic_rank(n, k, d, comp, ctx), rank);
    } else {
      rr = bic_rank(n, k, d, comp, ctx);
    }

    assert(r == rr);

    check_valid_bounded_composition(comp, n, k, d);
  }

  free(comp);
}

void bic_run_single_round_trip(const uint32_t n, const uint32_t k,
                               const uint32_t d, const bic_ctx_t ctx) {
  bic_run_round_trips(n, k, d, 1, NULL, NULL, NULL, NULL, ctx);
}

void bic_print_config(const uint32_t n, const uint32_t k, const uint32_t d,
                      const bic_ctx_t ctx) {
  const char *cache_name = bic_cache_names[ctx->cache_type];
  const char *order_name = bic_order_names[ctx->order];
  const char *unrank_alg_name = bic_unrank_alg_names[ctx->unrank_alg];
  const char *strategy_name = bic_strategy_names[ctx->strategy];
  printf("n=%5u, k=%5u, d=%5u, m=%5d, b=%5.0f, c=%5s, o=%5s, a=%7s, s=%6s, ", n,
         k, d, bits_fit_bic(n, k, d), BIT_LENGTH, cache_name, order_name,
         unrank_alg_name, strategy_name);
}

void bic_print_perf_stats(const uint32_t it, const long double utime,
                          const long double ucycles, const long double rtime,
                          const long double rcycles) {
  printf("i=%5u, unrank_avg_ns=%14.2Lf, unrank_avg_cyc=%14.2Lf, "
         "rank_avg_ns=%14.2Lf, rank_avg_cyc=%14.2Lf",
         it, utime / it, ucycles / it, rtime / it, rcycles / it);
}

uint8_t parse_uint16(void *value, const char *arg) {
  if (strchr(arg, '-')) {
    return 1;
  }
  *(uint16_t *)value = strtoul(arg, NULL, 0);
  return 0;
}

uint8_t parse_uint32(void *value, const char *arg) {
  if (strchr(arg, '-')) {
    return 1;
  }
  *(uint32_t *)value = strtoul(arg, NULL, 0);
  return 0;
}

uint8_t parse_string(void *value, const char *arg) {
  *(const char **)value = arg;
  return 0;
}

void print_help(const char *name, cli_option_t options[]) {
  printf("Usage: %s [OPTIONS]\n", name);
  for (size_t i = 0; options[i].def != NULL; ++i) {
    const cli_option_def_t *def = options[i].def;

    char flag[4];
    if (def->opt.val <= UINT8_MAX) {
      snprintf(flag, 4, "-%c,", def->opt.val);
    } else {
      snprintf(flag, 4, "   ");
    }

    char req[12];
    if (def->opt.has_arg == required_argument) {
      snprintf(req, 12, "=<%s>", def->type_description);
    }

    printf("  %s --%s%s\n        %s\n", flag, def->opt.name, req,
           def->description);

    if (def->num_sub_options == 0) {
      continue;
    }

    printf("        Available options for <%s> are:\n", def->type_description);
    for (size_t j = 0; j < def->num_sub_options; ++j) {
      printf("          * `%s` (%s)\n", def->sub_option_names[j],
             def->sub_option_infos[j]);
    }
  }
}

void parse_args(int32_t argc, char **argv, cli_option_t options[]) {
  size_t count = 0;
  while (options[count].def) {
    count++;
  }

  struct option *long_opt =
      (struct option *)malloc(sizeof(struct option) * (count + 1));
  char *short_opt = (char *)malloc(sizeof(char) * (count * 2 + 1));
  char *ptr = short_opt;

  for (size_t i = 0; i < count; ++i) {
    long_opt[i] = options[i].def->opt;
    if (!long_opt[i].val) {
      continue;
    }

    *ptr++ = long_opt[i].val;
    if (long_opt[i].has_arg == required_argument) {
      *ptr++ = ':';
    }
  }

  memset(&long_opt[count], 0, sizeof(struct option));
  *ptr = '\0';

  for (;;) {
    int c = getopt_long(argc, argv, short_opt, long_opt, NULL);
    if (c == -1) {
      break;
    }

    for (size_t i = 0; i < count; ++i) {
      if (c != long_opt[i].val) {
        continue;
      }

      const char *arg =
          long_opt[i].has_arg == required_argument ? optarg : NULL;
      if (options[i].def->parser(options[i].value, arg) != 0) {
        fprintf(stderr, "Invalid value for --%s: %s\n", long_opt[i].name, arg);
      }
      break;
    }
  }

  free(long_opt);
  free(short_opt);
}
