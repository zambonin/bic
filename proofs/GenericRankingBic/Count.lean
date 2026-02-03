/--
  Abstract combinatorial count function.
  Corresponds to the combinatorial count $\mathfrak{C}_{d}(n, k)$ in the paper.
  Represents the number of compositions of $n$ into $k$ parts (with part limits $d$ if applicable).
  This function serves as the basis for rank bounds.
-/
axiom count (n k d : Nat) : Nat

/--
  Base case for count: 0 parts.
  $\mathfrak{C}_{d}(n, 0) = 0$ for implicit $n > 0$.
-/
axiom count_zero : ∀ n d, count n 0 d = 0

/--
  Base case for count: 1 part.
  $\mathfrak{C}_{d}(n, 1) = 1$.
  There is exactly one composition of $n$ into 1 part: `[n]`.
-/
axiom count_one  : ∀ n d, count n 1 d = 1
