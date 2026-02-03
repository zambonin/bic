import GenericRankingBic.RankingCorrectness
import GenericRankingBic.StratAbstract

/-- Corresponds to $\mathsf{D}_{\mathsf{colex}}(k) = (k-1, 1)$. -/
def D_colex : Nat → Nat × Nat
  | 0 => (0, 0)
  | 1 => (1, 0)
  | Nat.succ (Nat.succ k) => (Nat.succ k, 1)

/-- Proof that $(k-1) + 1 = k$ for $k > 1$. -/
lemma D_colex_sum :
  ∀ k, k > 1 → (D_colex k).1 + (D_colex k).2 = k := by
  intro k hk
  cases k with
  | zero => cases Nat.not_lt_zero _ hk
  | succ k' =>
      cases k' with
      | zero => cases Nat.lt_asymm hk hk
      | succ k'' =>
          simp [D_colex]

/-- For $k > 1$, $k-1 < k$ and $1 < k$. -/
lemma D_colex_lt :
  ∀ k, k > 1 → (D_colex k).1 < k ∧ (D_colex k).2 < k := by
  intro k hk
  cases k with
  | zero => cases Nat.not_lt_zero _ hk
  | succ k' =>
      cases k' with
      | zero => cases Nat.lt_asymm hk hk
      | succ k'' =>
          simp [D_colex]

/-- Implementation of the standard colexicographical order. -/
noncomputable def Sigma_colex : Strategy :=
{
  decompose := D_colex
  decompose_sum := D_colex_sum
  decompose_lt  := D_colex_lt

  path := path_std
  path_nodup := path_std_nodup

  split := split_std
  combine := combine_std
  split_bound := split_std_bound

  combine_split := combine_split_std
  split_combine := split_combine_std
}

/-- Corollary of Theorem 3 applied to the $\Sigma_{\mathsf{colex}}$ strategy. -/
theorem rank_unrank_colex :
  ∀ n k d r,
    r < count n k d →
    rank Sigma_colex n k d (unrank Sigma_colex n k d r) = r :=
by
  intro n k d r hr
  exact rank_unrank_rec Sigma_colex n k d r hr
