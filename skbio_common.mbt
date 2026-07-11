///|
/// Shared utilities for scikit-bio modules.

///|
/// Sum of a Double array.
pub fn sum_double(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for x in arr {
    s = s + x
  }
  s
}

///|
/// Mean of a Double array (0.0 if empty).
pub fn mean_double(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    sum_double(arr) / arr.length().to_double()
  }
}

///|
/// Count elements satisfying a predicate.
pub fn count_if(arr : Array[Double], pred : (Double) -> Bool) -> Int {
  let mut n = 0
  for x in arr {
    if pred(x) {
      n = n + 1
    }
  }
  n
}

///|
/// Fisher-Yates in-place shuffle.
pub fn fisher_yates_shuffle(arr : Array[Int], rng : @random.Rand) -> Unit {
  let n = arr.length()
  for i in 0..<n {
    let j = i + rng.int(limit=n - i)
    let tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
  }
}

///|
/// Create a chacha8 RNG from an Int seed.
pub fn make_rng(seed : Int) -> @random.Rand {
  @random.Rand::chacha8(
    seed=Bytes::makei(32, fn(i) { ((seed >> (i % 4 * 8)) & 0xFF).to_byte() }),
  )
}

///|
/// Iterate over the upper triangle of a square matrix (i < j).
/// Calls f(i, j) for each pair.
pub fn upper_triangle_iter(n : Int, f : (Int, Int) -> Unit) -> Unit {
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      f(i, j)
      j = j + 1
    }
    i = i + 1
  }
}
