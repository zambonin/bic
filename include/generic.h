#ifndef GENERIC_H
#define GENERIC_H

#include "common.h"

typedef struct {
  uint32_t n;
  uint32_t k_l;
  uint32_t k_r;
  uint32_t d;
  uint32_t current_i;
  uint32_t val_1;
  uint32_t val_2;
} bic_path_state_t;

static inline void compute_bounds(uint32_t n, uint32_t k_l, uint32_t k_r,
                                  uint32_t d, uint32_t *min_wl,
                                  uint32_t *max_wl) {
  uint64_t max_wr = (uint64_t)k_r * d;
  *min_wl = ((uint64_t)n > max_wr) ? n - (uint32_t)max_wr : 0;

  uint64_t max_wl_cap = (uint64_t)k_l * d;
  *max_wl = ((uint64_t)n < max_wl_cap) ? n : (uint32_t)max_wl_cap;
}

void split_std(uintx rank, uintx vol_L, uintx vol_R, uint32_t w_L, uint32_t w_R,
               uintx *r_L, uintx *r_R);
uintx unsplit_std(uintx r_L, uintx r_R, uintx vol_L, uintx vol_R, uint32_t w_L,
                  uint32_t w_R);

typedef void (*bic_decomp_func_t)(uint32_t k, uint32_t *k_l, uint32_t *k_r);
typedef void (*bic_path_init_func_t)(bic_path_state_t *state, uint32_t n,
                                     uint32_t k_l, uint32_t k_r, uint32_t d);
typedef int (*bic_path_next_func_t)(bic_path_state_t *state, uint32_t *w_l);
typedef void (*bic_split_func_t)(uintx rank, uintx vol_L, uintx vol_R,
                                 uint32_t w_L, uint32_t w_R, uintx *r_L,
                                 uintx *r_R);
typedef uintx (*bic_unsplit_func_t)(uintx r_L, uintx r_R, uintx vol_L,
                                    uintx vol_R, uint32_t w_L, uint32_t w_R);

typedef struct {
  const char *name;
  bic_decomp_func_t decomp;
  bic_path_init_func_t path_init;
  bic_path_next_func_t path_next;
  bic_split_func_t split;
  bic_unsplit_func_t unsplit;
} bic_path_strategy_t;

void bic_generic_unrank(uint32_t *res, uint32_t n, uint32_t k, uint32_t d,
                        const uintx *rank, const bic_path_strategy_t *strat,
                        const bic_ctx_t ctx);

uintx bic_generic_rank(uint32_t n, uint32_t k, uint32_t d, const uint32_t *res,
                       const bic_path_strategy_t *strat, const bic_ctx_t ctx);

#endif
