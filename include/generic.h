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
} bic_iter_state_t;

static inline void compute_bounds(uint32_t n, uint32_t k_l, uint32_t k_r,
                                  uint32_t d, uint32_t *min_wl,
                                  uint32_t *max_wl) {

  uint32_t max_wr = (uint32_t)k_r * d;
  *min_wl = (n > max_wr) ? (n - max_wr) : 0;

  uint32_t max_wl_cap = (uint32_t)k_l * d;
  *max_wl = (n < max_wl_cap) ? n : (uint32_t)max_wl_cap;
}

void map_std(uintx rank, uintx vol_L, uintx vol_R, uint32_t w_L, uint32_t w_R,
             uintx *r_L, uintx *r_R);
uintx unmap_std(uintx r_L, uintx r_R, uintx vol_L, uintx vol_R, uint32_t w_L,
                uint32_t w_R);

typedef void (*bic_topo_func_t)(uint32_t k, uint32_t *k_l, uint32_t *k_r);
typedef void (*bic_iter_init_func_t)(bic_iter_state_t *state, uint32_t n,
                                     uint32_t k_l, uint32_t k_r, uint32_t d);
typedef int (*bic_iter_next_func_t)(bic_iter_state_t *state, uint32_t *w_l);
typedef void (*bic_map_func_t)(uintx rank, uintx vol_L, uintx vol_R,
                               uint32_t w_L, uint32_t w_R, uintx *r_L,
                               uintx *r_R);
typedef uintx (*bic_unmap_func_t)(uintx r_L, uintx r_R, uintx vol_L,
                                  uintx vol_R, uint32_t w_L, uint32_t w_R);

typedef struct {
  const char *name;
  bic_topo_func_t topo;
  bic_iter_init_func_t iter_init;
  bic_iter_next_func_t iter_next;
  bic_map_func_t map;
  bic_unmap_func_t unmap;
} bic_iter_strategy_t;

void bic_generic_unrank(uint32_t *res, uint32_t n, uint32_t k, uint32_t d,
                        const uintx *rank, const bic_iter_strategy_t *strat,
                        const bic_ctx_t ctx);

uintx bic_generic_rank(uint32_t n, uint32_t k, uint32_t d, const uint32_t *res,
                       const bic_iter_strategy_t *strat, const bic_ctx_t ctx);

#endif
