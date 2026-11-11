# =============================================================================
# [IMPROVEMENT-5] Fast exact scoring (Algorithm 1 of the paper) in O(T log T)
#
# The reference implementation `scoth_changes()` (file scoth_scoring_algorithm.R)
# needs, at EACH of its T-1 iterations,
#   * O(k) work over all k remaining cuts (compute_diff_costs, pmax, which.min),
#   * O(k) copies of two data.frames (rows are deleted with l$df[-(i+1), ]),
#   * 2 calls of the cost function, each costing O(segment length)
#     (and the linear cost even calls lm() + predict(), about 0.6-0.9 ms per
#     call: see results/cost_call_microbenchmark.csv).
# => total time at least O(T^2).
#
# This file computes EXACTLY the same scores with better data structures:
#   (1) the current cuts are kept in a doubly-linked list (vectors prv / nxt),
#       so removing a cut and finding its two neighbours is O(1);
#   (2) the scores are kept in a binary min-heap ordered by (score, index), so
#       "argmin" is O(1) and every update is O(log T).  Ties are broken by the
#       smallest index, exactly like which.min() in the original code;
#   (3) after the removal of a cut b with neighbours p < b < q, only the gains
#       of p and q change (G(c, S) only depends on the two segments around c),
#       so only 2 new segment costs are needed per iteration;
#   (4) segment costs: short segments (length <= exact_max_len) are computed
#       with a bit-identical copy of the original cost functions (no lm():
#       .lm.fit() gives the same floating-point numbers), long segments are
#       computed in O(1) from prefix sums (cumulative sums of y, y^2, t*y).
# Total: O(T log T) time and O(T) memory, pure R.
# =============================================================================

# -----------------------------------------------------------------------------
# [IMPROVEMENT-5] Names of the cost functions supported by the fast algorithm
# -----------------------------------------------------------------------------
fast_cost_names = function() {
  c("cusum_L2", "reg_cusum_L2", "cusum_L2_multivariate", "reg_cusum_L2_multivariate")
}

#' [IMPROVEMENT-5] Can the fast algorithm be used on this series?
#'
#' It needs at least 2 finite numeric values (a univariate vector, or a
#' data.frame whose columns are all numeric, as in the TCPD data). For any other
#' input `scoth_scoring_func` keeps the reference implementation.
#' @param x univariate (vector) or multivariate (data.frame) series
#' @return TRUE or FALSE
fast_input_ok = function(x) {
  if (is.data.frame(x)) {
    ncol(x) >= 1 && nrow(x) >= 2 && all(vapply(x, is.numeric, TRUE)) &&
      all(vapply(x, function(v) all(is.finite(v)), TRUE))
  } else {
    is.numeric(x) && is.null(dim(x)) && length(x) >= 2 && all(is.finite(x))
  }
}

#' [IMPROVEMENT-5] Find the name of a cost function given as a function object
#'
#' @param cost_func a cost function (closure) or a cost function name
#' @return one of fast_cost_names(), or NULL when the cost function is not
#' supported by the fast algorithm (e.g. the dummy c_box used in the tests)
fast_cost_name = function(cost_func) {
  if (is.character(cost_func)) {
    return(if (cost_func %in% fast_cost_names()) cost_func else NULL)
  }
  all_f = all_cost_functions()
  for (nm in fast_cost_names()) {
    if (identical(cost_func, all_f[[nm]], ignore.environment = TRUE)) return(nm)
  }
  NULL
}

# -----------------------------------------------------------------------------
# [IMPROVEMENT-5] Bit-identical copies ("replicas") of the original cost functions
#
# They perform the same floating-point operations in the same order as the
# original functions of cost_functions.R, hence return identical numbers
# (checked with identical() in tests/testthat/test-improvement5.R), but they
# avoid the expensive formula / model.frame / predict machinery of lm().
#   * lm(y ~ idx) calls lm.fit(X, y) with X = [1, idx], i.e. the C routine
#     Cdqrls (Householder QR); .lm.fit(X, y) calls the same routine.
#   * predict() of a univariate lm returns X %*% beta  (computed again here);
#     predict() of a multivariate lm ("mlm") returns fitted.values = Y - residuals.
# -----------------------------------------------------------------------------
replica_c_cusum_L2 = function(x, s, e) {
  x_cur = x[s:(e-1)]                      # same code as c_cusum_L2
  sum((x_cur - mean(x_cur))^2)
}

replica_c_reg_cusum_L2 = function(x, s, e) {
  if (e - s == 1) return(0)               # 1 point: lm gives a perfect fit, cost exactly 0
  idx = s:(e-1)
  x_cur = x[idx]
  X = cbind(1, as.double(idx))            # model.matrix(x_cur ~ idx)
  z = .lm.fit(X, as.double(x_cur))        # same QR as lm()
  if (z$rank < 2L) return(c_reg_cusum_L2(x, s, e))  # rare rank-deficient case: call the original
  sum((drop(X %*% z$coefficients) - x_cur)^2)       # predict(reg) = X %*% beta
}

replica_c_cusum_L2_multivariate = function(cols, s, e) {
  # `cols` = list of the columns of the data.frame, each with its ORIGINAL type:
  # mean() of an integer vector and of a double vector use different algorithms in R,
  # so an integer column (e.g. "Volume" of the apple dataset) must stay integer.
  if (e - s == 1) return(0)
  idx = s:(e-1)
  # same as sapply(features, function(f) x_cur[,f] - mean(x_cur[,f])) followed by
  # sum(.^2): the matrix is summed column after column, as this concatenation
  dev = unlist(lapply(cols, function(v) { v = v[idx]; v - mean(v) }), use.names = FALSE)
  sum(dev^2)
}

replica_c_reg_cusum_L2_multivariate = function(M, s, e) {
  if (e - s == 1) return(0)
  idx = s:(e-1)
  Y = M[idx, , drop = FALSE]
  X = cbind(1, as.double(idx))
  z = .lm.fit(X, Y)
  if (z$rank < 2L) return(c_reg_cusum_L2_multivariate(M, s, e))  # rare rank-deficient case: original
  sum(((Y - z$residuals) - Y)^2)          # predict.mlm = fitted.values = Y - residuals
}

# -----------------------------------------------------------------------------
# [IMPROVEMENT-5] Compensated prefix sums
#
# P = c(0, cumsum(v)) gives sum(v[s:(e-1)]) = P[e] - P[s] in O(1).
# Floating-point issue: P[e] and P[s] are rounded to ~1e-16 * |P| (the size of
# the WHOLE cumulative sum), which can be large compared with the sum over a
# short segment.  We therefore keep a second vector `lo` with the rounding
# error of each stored partial sum (like Kahan summation):
#     d_i = hi_i - hi_{i-1}  is the increment actually stored,
#     v_i - d_i              is the part of v_i that was lost by rounding,
#     lo  = cumsum(v - d)    accumulates these small lost parts.
# Then sum(v[s:(e-1)]) = (hi[e] - hi[s]) + (lo[e] - lo[s]) with a relative
# accuracy similar to a direct summation over the segment.
# -----------------------------------------------------------------------------
compensated_prefix_sums = function(v) {
  if (is.matrix(v)) {
    hi = rbind(0, apply(v, 2, cumsum))
    lo = rbind(0, apply(v - (hi[-1, , drop = FALSE] - hi[-nrow(hi), , drop = FALSE]), 2, cumsum))
  } else {
    hi = c(0, cumsum(v))
    lo = c(0, cumsum(v - diff(hi)))
  }
  list(hi = hi, lo = lo)
}

#' [IMPROVEMENT-5] Exact (bit-identical) segment cost function x[s:(e-1)] -> cost
#'
#' @param x univariate (vector) or multivariate (data.frame) series
#' @param cost_name one of fast_cost_names()
#' @return a function(s, e) returning the same number as the original cost function
fast_exact_cost = function(x, cost_name) {
  linear = grepl("^reg_", cost_name)
  if (grepl("_multivariate$", cost_name)) {
    if (linear) {
      M = as.matrix(x)                 # lm(as.matrix(x_cur) ~ idx) works on a double matrix
      storage.mode(M) = "double"
      function(s, e) replica_c_reg_cusum_L2_multivariate(M, s, e)
    } else {
      cols = lapply(seq_len(ncol(x)), function(f) x[[f]])   # columns with their own type
      function(s, e) replica_c_cusum_L2_multivariate(cols, s, e)
    }
  } else {
    if (linear) function(s, e) replica_c_reg_cusum_L2(x, s, e)
    else function(s, e) replica_c_cusum_L2(x, s, e)
  }
}

#' [IMPROVEMENT-5] Build the segment cost function used by the fast algorithm
#'
#' @param x univariate (vector) or multivariate (data.frame) series
#' @param cost_name one of fast_cost_names()
#' @param exact_max_len segments with at most this number of points use the
#' bit-identical replica of the original cost function (O(length) each); longer
#' segments use the O(1) prefix-sum formula. Inf = always exact, 0 = always O(1).
#' @return a list with `cost(s, e)` (hybrid cost of x[s:(e-1)]), `cost_exact(s, e)`,
#' `cost_prefix(s, e)` and `T` (length of the series)
fast_segment_cost = function(x, cost_name, exact_max_len = 64) {
  multivariate = grepl("_multivariate$", cost_name)
  linear = grepl("^reg_", cost_name)
  cost_exact = fast_exact_cost(x, cost_name)
  if (multivariate) {
    M = as.matrix(x)
    storage.mode(M) = "double"
  } else {
    M = matrix(as.double(x), ncol = 1)
  }
  T = nrow(M)

  # ---- O(1) costs from prefix sums --------------------------------------------
  # Centering (L2) or global detrending (linear cost) does not change the costs
  # mathematically (a constant, resp. a global line, is absorbed by the
  # segment mean, resp. the segment regression line), but it makes the numbers
  # small, which reduces cancellation in  sum(y^2) - sum(y)^2/n.
  # (Scaling would not help: floating-point relative errors are scale-free.)
  tt = seq_len(T) - (T + 1) / 2                          # centered time index
  if (linear) {
    Y = M - (cbind(1, tt) %*% .lm.fit(cbind(1, tt), M)$coefficients)  # global detrending
  } else {
    Y = M - rep(colMeans(M), each = T)                    # global centering
  }
  P1 = compensated_prefix_sums(Y)          # sum of y
  P2 = compensated_prefix_sums(Y * Y)      # sum of y^2
  P3 = if (linear) compensated_prefix_sums(Y * tt) else list()   # sum of t*y (linear cost only)
  H1 = P1$hi; L1 = P1$lo; H2 = P2$hi; L2 = P2$lo; H3 = P3$hi; L3 = P3$lo
  F = ncol(M)

  # constant runs: run_end[i] = last index of the run of identical rows containing i,
  # so that x[s:(e-1)] is constant  <=>  run_end[s] >= e-1  (cost exactly 0, as in the original)
  row_change = c(TRUE, if (T > 1) rowSums(M[-1, , drop = FALSE] != M[-T, , drop = FALSE]) > 0 else logical(0))
  run_id = cumsum(row_change)
  run_end = cumsum(rle(run_id)$lengths)[run_id]

  cost_prefix = function(s, e) {
    n = as.double(e - s)                   # double: n^3 would overflow the 32-bit integers
    if (n <= 1 || run_end[s] >= e - 1L) return(0)
    if (F == 1L) {
      S1 = (H1[e] - H1[s]) + (L1[e] - L1[s])
      S2 = (H2[e] - H2[s]) + (L2[e] - L2[s])
    } else {
      S1 = (H1[e, ] - H1[s, ]) + (L1[e, ] - L1[s, ])
      S2 = (H2[e, ] - H2[s, ]) + (L2[e, ] - L2[s, ])
    }
    Syy = S2 - S1 * S1 / n                 # sum of (y - mean)^2   (per feature)
    if (linear) {
      if (n == 2) return(0)                # a line through 2 points: perfect fit
      S3 = if (F == 1L) (H3[e] - H3[s]) + (L3[e] - L3[s]) else (H3[e, ] - H3[s, ]) + (L3[e, ] - L3[s, ])
      t_bar = (s + e - 1) / 2 - (T + 1) / 2      # mean of the centered times of the segment
      Sty = S3 - t_bar * S1                       # sum of (t - t_bar) * (y - mean)
      Stt = n * (n * n - 1) / 12                  # sum of (t - t_bar)^2 (consecutive integers)
      Syy = Syy - Sty * Sty / Stt                 # residual sum of squares of the line fit
    }
    sum(pmax(Syy, 0))                      # costs are >= 0; sum over features
  }

  L0 = exact_max_len
  cost = function(s, e) {
    if (e - s <= L0) cost_exact(s, e) else cost_prefix(s, e)
  }
  list(cost = cost, cost_exact = cost_exact, cost_prefix = cost_prefix, T = T)
}

#' [IMPROVEMENT-5] Fast exact version of `scoth_changes` (Algorithm 1)
#'
#' Same output as `scoth_changes(x, cost_func, s = 1, e = T + 1)`: the vector of
#' the T-1 scores (maximum normalized gains), but computed in O(T log T).
#' It is called by `scoth_scoring_func(..., algorithm = "fast")`, or by every
#' scoring of the pipeline after options(scoth.scoring_algorithm = "fast").
#'
#' @param x univariate (vector) or multivariate (data.frame) series
#' @param cost_name one of fast_cost_names(): "cusum_L2", "reg_cusum_L2",
#' "cusum_L2_multivariate" or "reg_cusum_L2_multivariate" (a name, not a function)
#' @param verbose whether to print the advance of the algorithm
#' @param exact_max_len segments with at most this number of points use the
#' bit-identical replicas of the original cost functions, longer segments the
#' O(1) prefix sums (see `fast_segment_cost`). Default: options(scoth.fast_exact_max_len),
#' or 64 when this option is not set; Inf gives scores bit-identical to `scoth_changes`.
#' @return vector of length T-1 with the score of each possible change point
#' @examples
#' set.seed(1)
#' x = c(rnorm(100), rnorm(100, 3))
#' score_fast = scoth_changes_fast(x, "cusum_L2")
#' score_original = scoth_scoring_func(x, "cusum_L2", verbose = FALSE)$changepoint_gain
#' max(abs(score_fast - score_original))   # 0, or a few 1e-16
#' which.max(score_fast)                    # 100: the change point
#' @export
scoth_changes_fast = function(x, cost_name, verbose = FALSE,
                              exact_max_len = getOption("scoth.fast_exact_max_len", 64)) {
  oracle = fast_segment_cost(x, cost_name, exact_max_len)
  seg_cost = oracle$cost
  T = oracle$T
  k = T - 1L                                   # number of possible cuts

  cost_initial = oracle$cost_exact(1L, T + 1L) # V(empty), same number as the original
  if (cost_initial == 0) {                     # constant series: no possible gain
    return(rep(0, k))
  }

  # ---- Positions --------------------------------------------------------------
  # A "boundary" is a time index b in 1..T+1: the segment [a, b) contains
  # x[a], ..., x[b-1]. Boundaries 1 and T+1 are the borders; b = 2..T are the
  # cuts. Cut b is the change point of index j = b - 1 in the output vector.
  # Doubly-linked list of the remaining boundaries:
  nxt = c(2:(T + 1L), NA_integer_)             # nxt[b] = next remaining boundary
  prv = c(NA_integer_, 1:T)                    # prv[b] = previous remaining boundary
  # segc[a]  = cost of the segment [a, nxt[a])            ("df" of the original)
  # pairc[b] = cost of the segment [prv[b], nxt[b])       ("df2" of the original)
  segc = numeric(T)                            # single points: cost exactly 0
  pairc = numeric(T)
  pairc[2:T] = vapply(2:T, function(b) seg_cost(b - 1L, b + 1L), 0)

  # ---- Scores: line 6 of Algorithm 1 at the first iteration --------------------
  # gain of cut b = V(S \ b) - V(S) = pairc[b] - segc[prv[b]] - segc[b]  (= pairc here)
  sc = pmax(((pairc[2:T] - segc[1:k]) - segc[2:T]) / cost_initial, rep(0, k))  # sc[j], j = b - 1

  # ---- Binary min-heap of the cuts, ordered by (score, index) -------------------
  # hp[1] is the cut with the smallest score (smallest index among ties, like which.min);
  # pos[j] = position of cut j in hp (0 once removed). A sorted array is a valid heap.
  hp = order(sc + 0, seq_len(k))               # (+0 turns a possible -0 into 0)
  pos = integer(k)
  pos[hp] = seq_len(k)
  hn = k                                       # current heap size

  todo = integer(3)                            # cuts whose heap position must be repaired
  newv = numeric(3)                            # their new scores
  for (loop in seq_len(k)) {
    if (verbose && k > 1000 && loop %% 1000 == 0) print(paste0(loop, "/", k))

    # -- line 7: remove the cut with the smallest score (root of the heap) --------
    j = hp[1L]
    pos[j] = 0L
    last = hp[hn]
    hn = hn - 1L
    nt = 0L
    if (hn > 0L) {                             # the last heap element goes to the root ...
      hp[1L] = last; pos[last] = 1L
      nt = 1L; todo[1L] = last; newv[1L] = sc[last]   # ... and is sifted down below
    }
    b = j + 1L                                 # boundary of the removed cut
    p = prv[b]; q = nxt[b]                     # its neighbours
    nxt[p] = q; prv[q] = p                     # unlink b
    segc[p] = pairc[b]                         # [p, b) U [b, q) = [p, q), already known

    # -- line 6: only the gains of the two neighbours change ----------------------
    if (p > 1L) {                              # p is a cut: new pair [prv[p], q)
      pp = prv[p]
      pc = seg_cost(pp, q)
      pairc[p] = pc
      g = ((pc - segc[pp]) - segc[p]) / cost_initial
      if (g > sc[p - 1L]) { nt = nt + 1L; todo[nt] = p - 1L; newv[nt] = g }
    }
    if (q <= T) {                              # q is a cut: new pair [p, nxt[q])
      qq = nxt[q]
      pc = seg_cost(p, qq)
      pairc[q] = pc
      g = ((pc - segc[p]) - segc[q]) / cost_initial
      if (g > sc[q - 1L]) { nt = nt + 1L; todo[nt] = q - 1L; newv[nt] = g }
    }

    # -- heap repair: the keys only INCREASE (max in line 6), so sift down ---------
    # at most 3 sift-downs per iteration: the last heap element moved to the root
    # and the 2 neighbours p-1, q-1 (only when their score increased)
    w = 0L
    while (w < nt) {
      w = w + 1L
      el = todo[w]
      se = newv[w]
      sc[el] = se
      i = pos[el]
      repeat {
        l = i + i
        if (l > hn) break
        ci = hp[l]; cs = sc[ci]
        r = l + 1L
        if (r <= hn) {
          ri = hp[r]; rs = sc[ri]
          if (rs < cs || (rs == cs && ri < ci)) { l = r; ci = ri; cs = rs }
        }
        if (cs < se || (cs == se && ci < el)) {
          hp[i] = ci; pos[ci] = i; i = l       # move the child up
        } else {
          break
        }
      }
      hp[i] = el; pos[el] = i
    }
  }
  sc
}
