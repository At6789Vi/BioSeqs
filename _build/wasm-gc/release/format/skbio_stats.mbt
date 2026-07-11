///|
/// scikit-bio Stats module - Statistical functions for bioinformatics.

///|
/// Randomly subsample from a vector of counts, with or without replacement.
/// Returns a subsampled vector where the sum equals n.
pub fn subsample_counts(
  counts : Array[Double],
  n : Int,
  replace : Bool,
  seed : Int,
) -> Array[Double] {
  let total = sum_double(counts)
  let total_int = total.to_int()
  if n < 0 {
    return Array::make(counts.length(), 0.0)
  }
  if !replace && n > total_int {
    return Array::make(counts.length(), 0.0)
  }
  let rng = make_rng(seed)
  let result = Array::make(counts.length(), 0.0)
  if replace {
    for _ in 0..<n {
      let r = rng.double() * total
      let mut cumulative = 0.0
      let mut idx = 0
      while idx < counts.length() && cumulative + counts[idx] < r {
        cumulative = cumulative + counts[idx]
        idx = idx + 1
      }
      if idx >= counts.length() {
        idx = counts.length() - 1
      }
      result[idx] = result[idx] + 1.0
    }
  } else {
    let pool : Array[Int] = Array::new()
    for i in 0..<counts.length() {
      let cnt = counts[i].to_int()
      for _ in 0..<cnt {
        pool.push(i)
      }
    }
    fisher_yates_shuffle(pool, rng)
    let take = if n < pool.length() { n } else { pool.length() }
    for k in 0..<take {
      let idx = pool[k]
      result[idx] = result[idx] + 1.0
    }
  }
  result
}

///|
/// Bonferroni correction. p_adj = p * m, clipped to [0, 1].
pub fn p_adjust_bonferroni(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length().to_double()
  let result = Array::make(pvalues.length(), 0.0)
  for i in 0..<pvalues.length() {
    let adj = pvalues[i] * m
    result[i] = if adj > 1.0 { 1.0 } else { adj }
  }
  result
}

///|
/// Stepwise p-value adjustment shared by Holm and Benjamini-Hochberg.
/// `factor_fn(i, n)` scales the sorted p-value at rank i; `descending` selects
/// iteration direction and monotonicity (Holm: forward/max, BH: backward/min).
fn stepwise_adjust(
  pvalues : Array[Double],
  factor_fn : (Int, Int) -> Double,
  descending : Bool,
  prev_init : Double,
) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return Array::new()
  }
  let indexed = Array::makei(n, fn(i) { (i, pvalues[i]) })
  sort_pairs_by_p(indexed, 0, n - 1)
  let adjusted_sorted = Array::make(n, 0.0)
  let mut prev = prev_init
  for k in 0..<n {
    let i = if descending { n - 1 - k } else { k }
    let mut adj = indexed[i].1 * factor_fn(i, n)
    if descending {
      if adj > prev {
        adj = prev
      }
    } else if adj < prev {
      adj = prev
    }
    if adj > 1.0 {
      adj = 1.0
    }
    adjusted_sorted[i] = adj
    prev = adj
  }
  let result = Array::make(n, 0.0)
  for k in 0..<n {
    result[indexed[k].0] = adjusted_sorted[k]
  }
  result
}

///|
/// Holm-Bonferroni (step-down) correction.
pub fn p_adjust_holm(pvalues : Array[Double]) -> Array[Double] {
  stepwise_adjust(pvalues, fn(i, n) { (n - i).to_double() }, false, 0.0)
}

///|
/// Benjamini-Hochberg (FDR) correction.
pub fn p_adjust_fdr_bh(pvalues : Array[Double]) -> Array[Double] {
  stepwise_adjust(
    pvalues,
    fn(i, n) { n.to_double() / (i + 1).to_double() },
    true,
    1.0,
  )
}

///|
/// Quick sort for (index, pvalue) pairs by p-value ascending.
fn sort_pairs_by_p(arr : Array[(Int, Double)], lo : Int, hi : Int) -> Unit {
  if lo < hi {
    let pivot = arr[hi].1
    let mut i = lo - 1
    for j in lo..<hi {
      if arr[j].1 <= pivot {
        i = i + 1
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
      }
    }
    let tmp = arr[i + 1]
    arr[i + 1] = arr[hi]
    arr[hi] = tmp
    let pi = i + 1
    sort_pairs_by_p(arr, lo, pi - 1)
    sort_pairs_by_p(arr, pi + 1, hi)
  }
}

///|
/// Result of a permutation test.
pub struct PermutationTestResult {
  statistic : Double
  p_value : Double
  permutations : Int
}

///|
/// Perform a permutation test on two independent samples.
/// p-value is proportion of permuted stats >= original stat (one-sided, greater).
pub fn permutation_test(
  sample1 : Array[Double],
  sample2 : Array[Double],
  test_fn : (Array[Double], Array[Double]) -> Double,
  permutations : Int,
  seed : Int,
) -> PermutationTestResult {
  let orig_stat = test_fn(sample1, sample2)
  if permutations <= 0 {
    return PermutationTestResult::{
      statistic: orig_stat,
      p_value: 0.0 / 0.0,
      permutations: 0,
    }
  }
  let n1 = sample1.length()
  let n2 = sample2.length()
  let n_total = n1 + n2
  let combined = Array::make(n_total, 0.0)
  for i in 0..<n1 {
    combined[i] = sample1[i]
  }
  for i in 0..<n2 {
    combined[n1 + i] = sample2[i]
  }
  let s1 = Array::make(n1, 0.0)
  let s2 = Array::make(n2, 0.0)
  let rng = make_rng(seed)
  let mut count_ge = 0
  for _ in 0..<permutations {
    // Fisher-Yates shuffle in-place on combined (Array[Double], kept inline).
    for j in 0..<n_total {
      let k = j + rng.int(limit=n_total - j)
      let tmp = combined[j]
      combined[j] = combined[k]
      combined[k] = tmp
    }
    for j in 0..<n1 {
      s1[j] = combined[j]
    }
    for j in 0..<n2 {
      s2[j] = combined[n1 + j]
    }
    let perm_stat = test_fn(s1, s2)
    if perm_stat >= orig_stat {
      count_ge = count_ge + 1
    }
  }
  let p_value = count_ge.to_double() / permutations.to_double()
  PermutationTestResult::{ statistic: orig_stat, p_value, permutations }
}

///|
/// Difference of means test statistic (mean1 - mean2).
pub fn diff_of_means(s1 : Array[Double], s2 : Array[Double]) -> Double {
  mean_double(s1) - mean_double(s2)
}

///|
/// A distance matrix (square, symmetric, zero diagonal).
pub struct DistanceMatrix {
  data : Array[Double]
  size : Int
  ids : Array[String]
}

///|
/// Create a DistanceMatrix from a flat array (row-major order).
pub fn DistanceMatrix::from_flat(
  data : Array[Double],
  size : Int,
  ids : Array[String],
) -> DistanceMatrix {
  DistanceMatrix::{ data, size, ids }
}

///|
/// Create a DistanceMatrix from a 2D array.
pub fn DistanceMatrix::from_2d(
  matrix : Array[Array[Double]],
  ids : Array[String],
) -> DistanceMatrix {
  let n = matrix.length()
  let data = Array::make(n * n, 0.0)
  for i in 0..<n {
    for j in 0..<n {
      data[i * n + j] = matrix[i][j]
    }
  }
  DistanceMatrix::{ data, size: n, ids }
}

///|
/// Get distance between sample i and sample j.
pub fn DistanceMatrix::get(self : DistanceMatrix, i : Int, j : Int) -> Double {
  self.data[i * self.size + j]
}

///|
/// Result of an ANOSIM test.
pub struct ANOSIMResult {
  r_statistic : Double
  p_value : Double
  permutations : Int
  num_groups : Int
  sample_size : Int
}

///|
/// Perform ANOSIM (Analysis of Similarities) test.
/// groups: array of group labels (integers), one per sample.
pub fn anosim(
  dm : DistanceMatrix,
  groups : Array[Int],
  permutations : Int,
  seed : Int,
) -> ANOSIMResult {
  let n = dm.size
  let num_groups = get_unique_groups(groups).length()
  let (distances, pair_indices) = anosim_precompute(dm)
  let ranks = rank_data(distances)
  let num_pairs = distances.length()
  let orig_r = anosim_r_from_ranks(ranks, pair_indices, groups, n, num_pairs)
  if permutations <= 0 {
    return ANOSIMResult::{
      r_statistic: orig_r,
      p_value: 0.0 / 0.0,
      permutations: 0,
      num_groups,
      sample_size: n,
    }
  }
  let rng = make_rng(seed)
  let mut count_ge = 0
  for _ in 0..<permutations {
    let perm_groups = groups.copy()
    fisher_yates_shuffle(perm_groups, rng)
    let perm_r = anosim_r_from_ranks(
      ranks, pair_indices, perm_groups, n, num_pairs,
    )
    if perm_r >= orig_r {
      count_ge = count_ge + 1
    }
  }
  let p_value = count_ge.to_double() / permutations.to_double()
  ANOSIMResult::{
    r_statistic: orig_r,
    p_value,
    permutations,
    num_groups,
    sample_size: n,
  }
}

///|
/// Precompute distances and pair index mapping (constant across permutations).
fn anosim_precompute(dm : DistanceMatrix) -> (Array[Double], Array[(Int, Int)]) {
  let n = dm.size
  let num_pairs = n * (n - 1) / 2
  let distances = Array::make(num_pairs, 0.0)
  let pair_indices = Array::make(num_pairs, (0, 0))
  let mut idx = 0
  upper_triangle_iter(n, fn(i, j) {
    distances[idx] = dm.get(i, j)
    pair_indices[idx] = (i, j)
    idx = idx + 1
  })
  (distances, pair_indices)
}

///|
/// Compute ANOSIM R statistic from precomputed ranks.
fn anosim_r_from_ranks(
  ranks : Array[Double],
  pair_indices : Array[(Int, Int)],
  groups : Array[Int],
  n : Int,
  num_pairs : Int,
) -> Double {
  let mut sum_within = 0.0
  let mut count_within = 0.0
  let mut sum_between = 0.0
  let mut count_between = 0.0
  for i in 0..<num_pairs {
    let (a, b) = pair_indices[i]
    if groups[a] != groups[b] {
      sum_between = sum_between + ranks[i]
      count_between = count_between + 1.0
    } else {
      sum_within = sum_within + ranks[i]
      count_within = count_within + 1.0
    }
  }
  let r_b = if count_between > 0.0 { sum_between / count_between } else { 0.0 }
  let r_w = if count_within > 0.0 { sum_within / count_within } else { 0.0 }
  let denom = n.to_double() * (n.to_double() - 1.0) / 4.0
  if denom == 0.0 {
    0.0
  } else {
    (r_b - r_w) / denom
  }
}

///|
/// Compute ranks of data (average ranks for ties).
fn rank_data(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let result = Array::make(n, 0.0)
  let indexed = Array::makei(n, fn(i) { (i, data[i]) })
  sort_pairs_by_p(indexed, 0, n - 1)
  let mut i = 0
  while i < n {
    let mut j = i
    while j + 1 < n && indexed[j + 1].1 == indexed[i].1 {
      j = j + 1
    }
    let avg_rank = (i.to_double() + j.to_double()) / 2.0 + 1.0
    for k in i..<(j + 1) {
      result[indexed[k].0] = avg_rank
    }
    i = j + 1
  }
  result
}

///|
/// Result of a PERMANOVA test.
pub struct PERMANOVAResult {
  f_statistic : Double
  p_value : Double
  permutations : Int
  num_groups : Int
  sample_size : Int
  df_between : Int
  df_within : Int
  ss_between : Double
  ss_within : Double
  ss_total : Double
}

///|
/// Perform PERMANOVA test using pseudo-F statistic on distance matrix.
pub fn permanova(
  dm : DistanceMatrix,
  groups : Array[Int],
  permutations : Int,
  seed : Int,
) -> PERMANOVAResult {
  let n = dm.size
  let ss_total = permanova_ss_total(dm)
  let group_ids = get_unique_groups(groups)
  let num_groups = group_ids.length()
  let (f_stat, ss_bet, ss_with, df_bet, df_with) = permanova_f_from_ss(
    dm, groups, group_ids, n, ss_total,
  )
  if permutations <= 0 {
    return PERMANOVAResult::{
      f_statistic: f_stat,
      p_value: 0.0 / 0.0,
      permutations: 0,
      num_groups,
      sample_size: n,
      df_between: df_bet,
      df_within: df_with,
      ss_between: ss_bet,
      ss_within: ss_with,
      ss_total,
    }
  }
  let rng = make_rng(seed)
  let mut count_ge = 0
  for _ in 0..<permutations {
    let perm_groups = groups.copy()
    fisher_yates_shuffle(perm_groups, rng)
    let (perm_f, _, _, _, _) = permanova_f_from_ss(
      dm, perm_groups, group_ids, n, ss_total,
    )
    if perm_f >= f_stat {
      count_ge = count_ge + 1
    }
  }
  let p_value = count_ge.to_double() / permutations.to_double()
  PERMANOVAResult::{
    f_statistic: f_stat,
    p_value,
    permutations,
    num_groups,
    sample_size: n,
    df_between: df_bet,
    df_within: df_with,
    ss_between: ss_bet,
    ss_within: ss_with,
    ss_total,
  }
}

///|
/// Precompute ss_total from distance matrix (constant across permutations).
fn permanova_ss_total(dm : DistanceMatrix) -> Double {
  let n = dm.size
  let mut ss_total = 0.0
  upper_triangle_iter(n, fn(i, j) {
    let d = dm.get(i, j)
    ss_total = ss_total + d * d
  })
  ss_total / n.to_double()
}

///|
/// Extract unique group IDs without HashMap (faster for small group counts).
fn get_unique_groups(groups : Array[Int]) -> Array[Int] {
  let group_ids : Array[Int] = Array::new()
  for g in groups {
    let mut found = false
    for existing in group_ids {
      if existing == g {
        found = true
        break
      }
    }
    if !found {
      group_ids.push(g)
    }
  }
  group_ids
}

///|
/// Compute PERMANOVA pseudo-F statistic using precomputed ss_total.
fn permanova_f_from_ss(
  dm : DistanceMatrix,
  groups : Array[Int],
  group_ids : Array[Int],
  n : Int,
  ss_total : Double,
) -> (Double, Double, Double, Int, Int) {
  let num_groups = group_ids.length()
  let mut ss_within = 0.0
  for g in group_ids {
    let mut ss_g = 0.0
    let mut n_g = 0
    for i in 0..<n {
      if groups[i] == g {
        for j in (i + 1)..<n {
          if groups[j] == g {
            let d = dm.get(i, j)
            ss_g = ss_g + d * d
          }
        }
        n_g = n_g + 1
      }
    }
    if n_g > 0 {
      ss_within = ss_within + ss_g / n_g.to_double()
    }
  }
  let ss_between = ss_total - ss_within
  let df_between = num_groups - 1
  let df_within = n - num_groups
  let f_stat = if df_within > 0 && ss_within > 0.0 {
    ss_between / df_between.to_double() / (ss_within / df_within.to_double())
  } else {
    0.0
  }
  (f_stat, ss_between, ss_within, df_between, df_within)
}

///|
/// Calculate a confidence bound assuming a normal distribution.
/// Returns (lower, upper) confidence interval bounds.
pub fn confidence_bound(
  mean : Double,
  std : Double,
  n : Int,
  confidence : Double,
) -> (Double, Double) {
  if n <= 1 || std == 0.0 {
    (mean, mean)
  } else {
    let z = normal_ppf(0.5 + confidence / 2.0)
    let se = std / n.to_double().sqrt()
    (mean - z * se, mean + z * se)
  }
}

///|
/// Evaluate a polynomial sum(coeffs[i] * x^i) via Horner's method.
fn poly_eval(coeffs : Array[Double], x : Double) -> Double {
  let mut r = 0.0
  for c in coeffs {
    r = r * x + c
  }
  r
}

///|
/// Approximate inverse normal CDF (probit function).
/// Uses Beasley-Springer-Moro algorithm.
fn normal_ppf(p : Double) -> Double {
  if p <= 0.0 {
    return -1000000000000000000000000000000.0
  }
  if p >= 1.0 {
    return 1000000000000000000000000000000.0
  }
  let a = [
    -39.69683028665376, 220.9460984245205, -275.9285104469687, 138.3577518672690,
    -30.66479806614716, 2.506628277459239,
  ]
  let b = [
    -54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972,
    -13.28068155288572,
  ]
  let c = [
    -0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734,
    4.374664141464968, 2.938163982698783,
  ]
  let d = [
    0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416,
  ]
  let p_low = 0.02425
  let p_high = 1.0 - p_low
  if p < p_low {
    let q = (-2.0 * @math.ln(p)).sqrt()
    poly_eval(c, q) / (poly_eval(d, q) * q + 1.0)
  } else if p <= p_high {
    let q = p - 0.5
    let r = q * q
    poly_eval(a, r) * q / (poly_eval(b, r) * r + 1.0)
  } else {
    let q = (-2.0 * @math.ln(1.0 - p)).sqrt()
    -poly_eval(c, q) / (poly_eval(d, q) * q + 1.0)
  }
}
