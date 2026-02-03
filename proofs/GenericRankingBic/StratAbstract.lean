import Mathlib.Tactic

set_option linter.unusedVariables false

/--
  Strategy $\Sigma = (D, P, S)$ for generic ranking.
  Defines the recursive decomposition structure required for ranking and unranking.

  Components:
  - `D` (Decomposition): Splits the problem size $k$ into $(k_L, k_R)$.
  - `P` (Path): Determines the sequence of block sizes for offset calculation.
  - `S` (Split): Splits the residual rank into $(r_L, r_R)$ based on subproblem capacities.

  The strategy must satisfy specific axioms (bounds, sums, and inverses) to guarantee correctness.
-/
structure Strategy where
  /-- Decomposition function $\mathsf{D}$. splits $k$ into $(k_L, k_R)$. -/
  decompose : Nat → Nat × Nat
  /-- The parts returned by decompose must sum to k. -/
  decompose_sum : ∀ k, k > 1 → (decompose k).1 + (decompose k).2 = k
  /-- The parts returned by decompose must be strictly smaller than k. -/
  decompose_lt  : ∀ k, k > 1 → (decompose k).1 < k ∧ (decompose k).2 < k

  /-- Path function $\mathsf{P}$. Determines the sequence of block sizes. -/
  path : Nat → Nat → Nat → List Nat
  /-- The path must not contain duplicates. -/
  path_nodup : ∀ n kl kr, (path n kl kr).Nodup

  /-- Split function $\mathsf{S}$. Splits a residual rank $r$ into $(r_L, r_R)$. -/
  split : (nl nr wl : Nat) → Nat → Nat × Nat
  /-- Combine function $\mathsf{S}^{-1}$. Combines $(r_L, r_R)$ into a residual rank. -/
  combine : (nl nr wl : Nat) → Nat → Nat → Nat
  /-- Ensures the split outputs respect the cardinalities of the subproblems. -/
  split_bound :
    ∀ nl nr wl r,
      r < nl * nr →
      (split nl nr wl r).1 < nl ∧ (split nl nr wl r).2 < nr
  /-- $\mathsf{combine}(\mathsf{split}(r)) = r$. -/
  combine_split : (nl nr wl : Nat) → ∀ r, combine nl nr wl (split nl nr wl r).1 (split nl nr wl r).2 = r
  /-- $\mathsf{split}(\mathsf{combine}(r_L, r_R)) = (r_L, r_R)$. -/
  split_combine : (nl nr wl : Nat) → ∀ rl rr, rl < nl → rr < nr → split nl nr wl (combine nl nr wl rl rr) = (rl, rr)

/-- Corresponds to the increasing path $\mathsf{P}_{\mathsf{inc}}$. -/
noncomputable def path_std (n kl kr : Nat) : List Nat :=
  List.range (n + 1)

/-- Proof that path_std produces distinct elements (range is always distinct). -/
lemma path_std_nodup :
  ∀ n kl kr, (path_std n kl kr).Nodup := by
  intros; unfold path_std; exact List.nodup_range

/-- Default split function $\mathsf{S}(r)$. -/
def split_std (nl nr wl r : Nat) : Nat × Nat :=
  (r / nr, r % nr)

/-- Default inverse split function $\mathsf{S}^{-1}(r_L, r_R)$. -/
def combine_std (nl nr wl rl rr : Nat) : Nat :=
  rl * nr + rr

/-- Proof that combine_std is the left inverse of split_std. -/
lemma combine_split_std (nl nr wl : Nat) :
  ∀ r, combine_std nl nr wl (split_std nl nr wl r).1 (split_std nl nr wl r).2 = r := by
  intro r
  unfold split_std combine_std
  dsimp
  rw [Nat.mul_comm, Nat.div_add_mod]

/-- Proof that split_std is the left inverse of combine_std. -/
lemma split_combine_std (nl nr wl : Nat) :
  ∀ rl rr,
    rl < nl →
    rr < nr →
    split_std nl nr wl (combine_std nl nr wl rl rr) = (rl, rr) := by
  intro rl rr hrl hrr
  unfold split_std combine_std
  have hnr : nr > 0 := Nat.zero_lt_of_lt hrr
  apply Prod.ext
  · dsimp
    rw [Nat.add_comm, Nat.add_mul_div_right _ _ hnr, Nat.div_eq_of_lt hrr, Nat.zero_add]
  · dsimp
    rw [Nat.mul_add_mod_self_right, Nat.mod_eq_of_lt hrr]

/-- If r < nl * nr, then (r / nr) < nl and (r % nr) < nr. -/
lemma split_std_bound (nl nr wl : Nat) :
  ∀ r, r < nl * nr → (split_std nl nr wl r).1 < nl ∧ (split_std nl nr wl r).2 < nr := by
  intro r hr
  unfold split_std
  constructor
  · -- r / nr < nl
    rw [Nat.mul_comm] at hr
    exact Nat.div_lt_of_lt_mul hr
  · -- r % nr < nr
    have h_nl_nr_pos : 0 < nl * nr := Nat.zero_lt_of_lt hr
    have h_nr_pos : 0 < nr := by
      cases nr with
      | zero => rw [Nat.mul_zero] at h_nl_nr_pos; cases Nat.lt_irrefl 0 h_nl_nr_pos
      | succ n => exact Nat.zero_lt_succ n
    exact Nat.mod_lt r h_nr_pos
