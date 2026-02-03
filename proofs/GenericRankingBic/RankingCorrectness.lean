import GenericRankingBic.Count
import GenericRankingBic.RankingFunctions
import GenericRankingBic.StratAbstract

/-- Theorem 3: For all valid $r$, $\rank{\Sigma}(\unrank{\Sigma}(r)) = r$. -/
theorem rank_unrank_rec :
  ∀ (σ : Strategy) n k d r,
    r < count n k d →
    rank σ n k d (unrank σ n k d r) = r
:= by
  intro σ n k d r hr
  revert n d r
  induction k using Nat.strongRecOn with
  | ind k ih =>
      intro n d r hr
      cases k with
      | zero =>
          simp [rank, unrank]
          rw [count_zero] at hr
          cases Nat.not_lt_zero _ hr
      | succ k' =>
          cases k' with
          | zero =>
              simp [rank, unrank]
              rw [count_one] at hr
              have : r = 0 := Nat.eq_zero_of_le_zero (Nat.le_of_lt_succ hr)
              subst this
              rfl
          | succ k'' =>
              let k_curr := k'' + 2
              have h_k_gt_1 : k_curr > 1 := by omega

              let kl := (σ.decompose k_curr).1
              let kr := (σ.decompose k_curr).2
              let path := σ.path n kl kr
              let fb := find_block n kl kr d path r
              let wl := fb.1
              let r' := fb.2

              have hfb : r = offset n kl kr d path wl + r' :=
                find_block_correct n kl kr d path r wl r' rfl

              have hfb_bound : r' < count wl kl d * count (n - wl) kr d :=
                find_block_bound n kl kr d path r wl r' rfl

              let s := σ.split (count wl kl d) (count (n - wl) kr d) wl r'

              have hs_bound : s.1 < count wl kl d ∧ s.2 < count (n - wl) kr d :=
                Strategy.split_bound σ (count wl kl d) (count (n - wl) kr d) wl r' hfb_bound

              have hl : rank σ wl kl d (unrank σ wl kl d s.1) = s.1 :=
                ih kl (σ.decompose_lt _ h_k_gt_1).1 wl d s.1 hs_bound.1

              have hr_proof : rank σ (n - wl) kr d (unrank σ (n - wl) kr d s.2) = s.2 :=
                ih kr (σ.decompose_lt _ h_k_gt_1).2 (n - wl) d s.2 hs_bound.2

              have husum : wl = (unrank σ wl kl d s.1).foldl Nat.add 0 :=
                (unrank_sum σ wl kl d s.1).symm

              have h_sum : kl + kr = k_curr := σ.decompose_sum k_curr h_k_gt_1

              -- Expansion and alignment
              unfold unrank
              dsimp only
              have h_k_align : (k'' + 1 + 1) = k_curr := rfl
              rw [h_k_align, ← h_sum]

              -- Apply rank_decompose
              rw [rank_decompose σ n kl kr d wl]
              · -- Final identity
                rw [hl, hr_proof, σ.combine_split]
                exact hfb.symm
              · exact husum
