#include <assert.h>
#include <errno.h>

#include "extra.h"

static bool fail_malloc = false;
static bool fail_calloc = false;

#ifdef __cplusplus
extern "C" {
#endif
void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);

void *__wrap_malloc(size_t size) {
  if (fail_malloc) {
    fail_malloc = false;
    return NULL;
  }
  return __real_malloc(size);
}

void *__wrap_calloc(size_t count, size_t size) {
  if (fail_calloc) {
    fail_calloc = false;
    return NULL;
  }
  return __real_calloc(count, size);
}
#ifdef __cplusplus
}
#endif

static uintx wrong_rank(const uint32_t n, const uint32_t k, const uint32_t d,
                        const uint32_t *parts, const bic_ctx_t ctx) {
  (void)n;
  (void)k;
  (void)d;
  (void)parts;
  (void)ctx;
  return 1;
}

int32_t main(void) {
  bic_ctx_t ctx = bic_ctx_init();
  assert(ctx != NULL);
  assert(bic_ctx_set_order_by_name(NULL, ctx) == BIC_ERROR_ARGUMENT);
  assert(bic_ctx_set_order_by_name("missing", ctx) == BIC_ERROR_ARGUMENT);
  assert(bic_ctx_set_order(BIC_ORDER_COLEX, NULL) == BIC_ERROR_ARGUMENT);
  assert(bic_ctx_set_defaults(NULL) == BIC_ERROR_ARGUMENT);
  assert(bic_precompute(0, 2, 1, NULL) == BIC_ERROR_ARGUMENT);
  assert(bic_precompute(0, 2, 1, ctx) == BIC_SUCCESS);
  assert(bic_precompute(UINT32_MAX, 2, 1, ctx) == BIC_ERROR_RANGE);
  assert(errno == ERANGE && ctx->prepared && ctx->precomputed_n == 0);

  fail_malloc = true;
  long double time = 17;
  assert(bic_run_round_trips(0, 2, 1, 1, &time, NULL, NULL, NULL, ctx) ==
         BIC_ERROR_ALLOCATION);
  assert(!fail_malloc && time == 17);
  ctx->rank = wrong_rank;
  assert(bic_run_round_trips(0, 2, 1, 1, NULL, NULL, NULL, NULL, ctx) ==
         BIC_ERROR_VALIDATION);

  assert(bic_ctx_set_defaults(ctx) == BIC_SUCCESS);
  fail_calloc = true;
  assert(bic_precompute(4, 3, 2, ctx) == BIC_ERROR_ALLOCATION);
  assert(!fail_calloc && errno == ENOMEM && !ctx->prepared);
  assert(ctx->bin_cache == NULL && ctx->comb_cache == NULL);
  assert(bic_precompute(4, 3, 2, ctx) == BIC_SUCCESS);
  bic_ctx_destroy(ctx);

  uint32_t value = 0;
  cli_option_t options[] = {{&opt_sum, &value}, {NULL, NULL}};
  char program[] = "errors";
  char bad[] = "-nno";
  char good[] = "-n2";
  char *args[] = {program, bad, good, NULL};
  optind = 0;
  opterr = 1;
  assert(parse_args(3, args, options) == BIC_ERROR_ARGUMENT);
  assert(value == 0 && opterr == 1);
  optind = 0;
  fail_malloc = true;
  assert(parse_args(1, args, options) == BIC_ERROR_ALLOCATION);
  assert(!fail_malloc);
  assert(parse_uint32(NULL, "0") == BIC_ERROR_ARGUMENT);
  assert(parse_string(&value, NULL) == BIC_ERROR_ARGUMENT);
  return BIC_SUCCESS;
}
