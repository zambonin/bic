#include <stdio.h>

#include "extra.h"

int32_t main(int32_t argc, char **argv) {
  bool error = false;
  uint16_t n = 0;
  uint16_t k = 0;
  uint16_t d = 0;
  char *order = NULL;

  cli_option_t options[] = {
      {&opt_sum, &n},       {&opt_parts, &k}, {&opt_bound, &d},
      {&opt_order, &order}, {NULL, NULL},
  };

  parse_args(argc, argv, options);

  bic_ctx_t ctx = bic_ctx_init();
  uint32_t *comp = (uint32_t *)malloc(k * sizeof(uint32_t));

  if (k == 0 || order == NULL || (uint64_t)n > (uint64_t)k * d || ctx == NULL ||
      bic_ctx_set_order_by_name(order, ctx) || comp == NULL) {
    print_help(argv[0], options);
    error = true;
    goto cleanup;
  }

  for (uintx r = 0; r < ctx->comp(n, k, d, ctx); ++r) {
    bic_unrank(comp, n, k, d, r, ctx);
    for (uint16_t i = 0; i < k; ++i) {
      printf("%u ", comp[i]);
      comp[i] = 0;
    }
    printf("\n");
  }

cleanup:
  free(comp);
  bic_ctx_destroy(ctx);

  return error;
}
