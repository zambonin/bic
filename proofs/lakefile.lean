import Lake
open Lake DSL

package "GenericRankingBic" where

require "leanprover-community" / "mathlib"

@[default_target]
lean_lib "GenericRankingBic" where
  roots := #[`GenericRankingBic]
