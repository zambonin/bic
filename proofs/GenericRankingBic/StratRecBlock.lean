import GenericRankingBic.RankingCorrectness
import GenericRankingBic.StratAbstract

/-- Corresponds to $\mathsf{D}_{\mathsf{rb}}(k) = (\lfloor k/2 \rfloor, \lceil k/2 \rceil)$. -/
def D_split : Nat → Nat × Nat
  | 1 => (1, 0)
  | k => (k / 2, k - k / 2)

/-- Proof that $\lfloor k/2 \rfloor + \lceil k/2 \rceil = k$ for any $k$. -/
lemma D_split_sum :
  ∀ k, k > 1 → (D_split k).1 + (D_split k).2 = k := by
  intro k hk
  unfold D_split
  cases k with
  | zero => omega
  | succ k' =>
      cases k' with
      | zero => omega
      | succ k'' =>
          simp
          omega

/-- Proof that `D_split` parts are strictly less than `k` for $k > 1$. -/
lemma D_split_lt :
  ∀ k, k > 1 → (D_split k).1 < k ∧ (D_split k).2 < k := by
  intro k hk
  unfold D_split
  cases k with
  | zero => omega
  | succ k' =>
      cases k' with
      | zero => omega
      | succ k'' =>
          simp
          omega

/-- Implementation of the recursive block order. -/
noncomputable def Sigma_rb : Strategy :=
{
  decompose := D_split
  decompose_sum := D_split_sum
  decompose_lt  := D_split_lt

  path := path_std
  path_nodup := path_std_nodup

  split := split_std
  combine := combine_std
  split_bound := split_std_bound

  combine_split := combine_split_std
  split_combine := split_combine_std
}

/-- Corollary of Theorem 3 applied to the $\Sigma_{\mathsf{rb}}$ strategy. -/
theorem rank_unrank_rb :
  ∀ n k d r,
    r < count n k d →
    rank Sigma_rb n k d (unrank Sigma_rb n k d r) = r :=
by
  intro n k d r hr
  exact rank_unrank_rec Sigma_rb n k d r hr
