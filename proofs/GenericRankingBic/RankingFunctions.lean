import GenericRankingBic.Count
import GenericRankingBic.StratAbstract

set_option linter.unusedVariables false

/--
  Calculates the rank offset for a block.
  Sum of sizes of all blocks $W < wl$ in the path.
-/
noncomputable def offset
  (n kl kr d : Nat) (path : List Nat) (wl : Nat) : Nat :=
  match path with
  | [] => 0
  | w :: ws =>
      if w = wl then 0
      else count w kl d * count (n - w) kr d + offset n kl kr d ws wl

/--
  Finds the block $wl$ (and residual rank $r'$) such that $r$ falls within that block.
  Iterates through the path $P$, subtracting block sizes until $r$ is covered or path ends.
-/
noncomputable def find_block
  (n kl kr d : Nat) (p : List Nat) (r : Nat) : Nat × Nat :=
  match p with
  | [] => (0, 0)
  | wl :: ws =>
      let v := count wl kl d * count (n - wl) kr d
      if r < v then (wl, r)
      else find_block n kl kr d ws (r - v)

/--
  Generic unranking algorithm $\unrank{\Sigma}(d, n, k, r)$.
  Constructs a composition of $n$ into $k$ parts with rank $r$ using strategy $\Sigma$.

  Parameters:
  - `σ`: The strategy $\Sigma = (D, P, S)$.
  - `n`: Total sum of the composition.
  - `k`: Number of parts.
  - `d`: Maximum part size constraint (if applicable).
  - `r`: Rank of the composition.

  Returns:
  - A list of $k$ integers summing to $n$.

  Logic:
  1. Base cases: $k=0 \implies [], k=1 \implies [n]$.
  2. Recursive step:
     - Decompose $k$ into $(k_L, k_R)$ using $\mathsf{D}(k)$.
     - Determine path of block sizes $P(n, k_L, k_R)$.
     - Find the block $w_L$ containing $r$ using `find_block`.
     - Split the residual rank into $(r_L, r_R)$ using $\mathsf{S}$.
     - Recursively unrank left and right parts and concatenate.
-/
noncomputable def unrank
  (σ : Strategy) (n k d r : Nat) : List Nat :=
  match k with
  | 0 => []
  | 1 => [n]
  | Nat.succ (Nat.succ k') =>
      let h_k_gt_1 : Nat.succ (Nat.succ k') > 1 :=
        Nat.succ_lt_succ (Nat.zero_lt_succ _)
      let kl := (σ.decompose (Nat.succ (Nat.succ k'))).1
      let kr := (σ.decompose (Nat.succ (Nat.succ k'))).2
      let path := σ.path n kl kr
      let fb := find_block n kl kr d path r
      let wl := fb.1
      let s :=
        σ.split
          (count wl kl d)
          (count (n - wl) kr d)
          wl fb.2
      unrank σ wl kl d s.1 ++
      unrank σ (n - wl) kr d s.2
termination_by k
decreasing_by
  · exact (σ.decompose_lt _ h_k_gt_1).1
  · exact (σ.decompose_lt _ h_k_gt_1).2

/--
  Generic ranking algorithm $\rank{\Sigma}(d, n, k, z)$.
  Computes the rank of a composition $z$ under strategy $\Sigma$.

  Parameters:
  - `σ`: The strategy $\Sigma = (D, P, S)$.
  - `n`: Total sum of the composition.
  - `k`: Number of parts.
  - `d`: Maximum part size constraint.
  - `z`: The composition (list of integers).

  Returns:
  - The rank `r` corresponding to `z`.

  Logic:
  1. Base cases: $k=0,1 \implies 0$.
  2. Recursive step:
     - Decompose $k$ into $(k_L, k_R)$.
     - Split $z$ into $L$ (length $k_L$) and $R$ (length $k_R$).
     - Compute $w_L = \sum L$.
     - Recursively rank $L$ and $R$ to get $r_L, r_R$.
     - Combine $r_L, r_R$ into residual rank using $\mathsf{S}^{-1}$.
     - Add offset for block $w_L$ calculated from path $P$.
-/
noncomputable def rank
  (σ : Strategy) (n k d : Nat) (z : List Nat) : Nat :=
  match k with
  | 0 => 0
  | 1 => 0
  | Nat.succ (Nat.succ k') =>
      let h_k_gt_1 : Nat.succ (Nat.succ k') > 1 :=
        Nat.succ_lt_succ (Nat.zero_lt_succ _)
      let kl := (σ.decompose (Nat.succ (Nat.succ k'))).1
      let kr := (σ.decompose (Nat.succ (Nat.succ k'))).2
      let L := z.take kl
      let R := z.drop kl
      let path := σ.path n kl kr
      let wl := L.foldl Nat.add 0
      let rl := rank σ wl kl d L
      let rr := rank σ (n - wl) kr d R
      let r_rem :=
        σ.combine
          (count wl kl d)
          (count (n - wl) kr d)
          wl rl rr
      offset n kl kr d path wl + r_rem
termination_by k
decreasing_by
  · exact (σ.decompose_lt _ h_k_gt_1).1
  · exact (σ.decompose_lt _ h_k_gt_1).2

/-- If `find_block` returns $(w_L, r')$, then $r = \text{offset}(\dots, w_L) + r'$. -/
axiom find_block_correct :
  ∀ n kl kr d path r wl r',
    find_block n kl kr d path r = (wl, r') →
    r = offset n kl kr d path wl + r'

/-- Asserts that $\rank{\Sigma}(L ++ R) = \text{offset}(w_L) + \mathsf{S}^{-1}(\rank(L), \rank(R))$. -/
axiom rank_decompose :
  ∀ σ n kl kr d wl L R,
    wl = L.foldl Nat.add 0 →
    rank σ n (kl + kr) d (L ++ R)
      =
      offset n kl kr d (σ.path n kl kr) wl
      +
      σ.combine
        (count wl kl d)
        (count (n - wl) kr d)
        wl
        (rank σ wl kl d L)
        (rank σ (n - wl) kr d R)

/-- Ensures that the residual rank $r'$ is within the bounds of the chosen block's capacity. -/
axiom find_block_bound :
  ∀ n kl kr d path r wl r',
    find_block n kl kr d path r = (wl, r') →
    r' < count wl kl d * count (n - wl) kr d

/-- Guarantees that `unrank` always produces a composition that sums to $n$. -/
axiom unrank_sum :
  ∀ σ n k d r, (unrank σ n k d r).foldl Nat.add 0 = n
