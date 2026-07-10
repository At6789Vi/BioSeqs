///|
/// Error raised by skbio alignment operations.
pub suberror AlignError {
  AlignError(String)
}

///|
/// DNA sequence type.
pub struct DNA {
  seq : Seq
} derive(Eq)

///|
/// RNA sequence type.
pub struct RNA {
  seq : Seq
} derive(Eq)

///|
/// Protein sequence type.
pub struct Protein {
  seq : Seq
} derive(Eq)

///|
/// Create a DNA sequence from a string.
pub fn DNA::new(s : String) -> DNA {
  DNA::{ seq: Seq::new(s) }
}

///|
/// Create an RNA sequence from a string.
pub fn RNA::new(s : String) -> RNA {
  RNA::{ seq: Seq::new(s) }
}

///|
/// Create a Protein sequence from a string.
pub fn Protein::new(s : String) -> Protein {
  Protein::{ seq: Seq::new(s) }
}

///|
/// Get the underlying sequence.
pub fn DNA::to_seq(self : DNA) -> Seq {
  self.seq
}

///|
/// Get the underlying sequence.
pub fn RNA::to_seq(self : RNA) -> Seq {
  self.seq
}

///|
/// Get the underlying sequence.
pub fn Protein::to_seq(self : Protein) -> Seq {
  self.seq
}

///|
/// Convert to string.
pub fn DNA::to_string(self : DNA) -> String {
  self.seq.to_string()
}

///|
/// Convert to string.
pub fn RNA::to_string(self : RNA) -> String {
  self.seq.to_string()
}

///|
/// Convert to string.
pub fn Protein::to_string(self : Protein) -> String {
  self.seq.to_string()
}

///|
/// Get length.
pub fn DNA::length(self : DNA) -> Int {
  self.seq.length()
}

///|
/// Get length.
pub fn RNA::length(self : RNA) -> Int {
  self.seq.length()
}

///|
/// Get length.
pub fn Protein::length(self : Protein) -> Int {
  self.seq.length()
}

///|
/// Get character at position.
pub fn DNA::get_char(self : DNA, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// Get character at position.
pub fn RNA::get_char(self : RNA, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// Get character at position.
pub fn Protein::get_char(self : Protein, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// TabularMSA stores a multiple sequence alignment in tabular form.
pub struct TabularMSA {
  records : Array[SeqRecord]
  dtype : String
} derive(Eq)

///|
/// Create a new TabularMSA from SeqRecords.
pub fn TabularMSA::new(
  records : Array[SeqRecord],
  dtype? : String = "DNA",
) -> TabularMSA raise AlignError {
  if records.length() > 0 {
    let first_len = records[0].length()
    for rec in records {
      if rec.length() != first_len {
        raise AlignError("All sequences must have the same length in an MSA")
      }
    }
  }
  TabularMSA::{ records, dtype }
}

///|
/// Create TabularMSA from DNA sequences.
pub fn TabularMSA::from_dna(sequences : Array[DNA]) -> TabularMSA {
  let records = Array::new()
  let mut i = 0
  for dna in sequences {
    records.push(SeqRecord::new(dna.to_seq(), id="seq" + i.to_string()))
    i = i + 1
  }
  TabularMSA::{ records, dtype: "DNA" }
}

///|
/// Create TabularMSA from RNA sequences.
pub fn TabularMSA::from_rna(sequences : Array[RNA]) -> TabularMSA {
  let records = Array::new()
  let mut i = 0
  for rna in sequences {
    records.push(SeqRecord::new(rna.to_seq(), id="seq" + i.to_string()))
    i = i + 1
  }
  TabularMSA::{ records, dtype: "RNA" }
}

///|
/// Create TabularMSA from Protein sequences.
pub fn TabularMSA::from_protein(sequences : Array[Protein]) -> TabularMSA {
  let records = Array::new()
  let mut i = 0
  for prot in sequences {
    records.push(SeqRecord::new(prot.to_seq(), id="seq" + i.to_string()))
    i = i + 1
  }
  TabularMSA::{ records, dtype: "Protein" }
}

///|
/// Get the number of sequences (rows).
pub fn TabularMSA::num_records(self : TabularMSA) -> Int {
  self.records.length()
}

///|
/// Get the number of columns.
pub fn TabularMSA::get_alignment_length(self : TabularMSA) -> Int {
  if self.records.length() == 0 {
    0
  } else {
    self.records[0].length()
  }
}

///|
/// Get the shape (rows, columns).
pub fn TabularMSA::shape(self : TabularMSA) -> (Int, Int) {
  (self.num_records(), self.get_alignment_length())
}

///|
/// Access the SeqRecord at index.
pub fn TabularMSA::get(self : TabularMSA, i : Int) -> SeqRecord? {
  if i >= 0 && i < self.records.length() {
    Some(self.records[i])
  } else {
    None
  }
}

///|
/// Get the column at position as a string.
pub fn TabularMSA::column(self : TabularMSA, col : Int) -> String? {
  if self.records.length() == 0 {
    return None
  }
  if col < 0 || col >= self.get_alignment_length() {
    return None
  }
  let buf = StringBuilder::new(size_hint=self.records.length())
  for rec in self.records {
    match rec.seq.get_char(col) {
      Some(c) => buf.write_char(c)
      None => buf.write_char('?')
    }
  }
  Some(buf.to_string())
}

///|
/// Compute consensus sequence.
pub fn TabularMSA::consensus(self : TabularMSA) -> Seq {
  let n_cols = self.get_alignment_length()
  let n_rows = self.num_records()
  if n_cols == 0 || n_rows == 0 {
    return Seq::new("")
  }

  let buf = StringBuilder::new(size_hint=n_cols)
  for col = 0; col < n_cols; col = col + 1 {
    let counts = Map([], capacity=4)
    for rec in self.records {
      match rec.seq.get_char(col) {
        Some(c) => {
          let key = c.to_string()
          counts[key] = counts.get(key).map(fn(v) { v + 1 }).unwrap_or(1)
        }
        None => ()
      }
    }

    let mut max_count = 0
    let mut consensus_char = '-'
    for k, v in counts {
      if v > max_count || (v == max_count && k > consensus_char.to_string()) {
        max_count = v
        consensus_char = k.unsafe_get(0).unsafe_to_char()
      }
    }
    buf.write_char(consensus_char)
  }
  Seq::new(buf.to_string())
}

///|
/// Compute conservation score for each column.
pub fn TabularMSA::conservation(self : TabularMSA) -> Array[Double] {
  let n_cols = self.get_alignment_length()
  let n_rows = self.num_records()
  if n_cols == 0 || n_rows == 0 {
    return Array::new()
  }

  let result = Array::new()
  for col = 0; col < n_cols; col = col + 1 {
    let mut count_a = 0
    let mut count_c = 0
    let mut count_g = 0
    let mut count_t = 0

    for rec in self.records {
      let s = rec.seq.to_string()
      if col < s.length() {
        let c = s.unsafe_get(col).unsafe_to_char()
        match c {
          'A' | 'a' => count_a = count_a + 1
          'C' | 'c' => count_c = count_c + 1
          'G' | 'g' => count_g = count_g + 1
          'T' | 't' | 'U' | 'u' => count_t = count_t + 1
          _ => ()
        }
      }
    }

    let max_count = count_a.max(count_c).max(count_g).max(count_t)
    let max_freq = max_count.to_double() / n_rows.to_double()
    result.push(max_freq)
  }
  result
}

///|
/// Alignment path operation type.
pub enum AlignOp {
  Match
  Insert
  Delete
  Gap
}

///|
/// Alignment path representing the operations to transform one sequence to another.
pub struct AlignPath {
  operations : Array[AlignOp]
  score : Double
  start1 : Int
  end1 : Int
  start2 : Int
  end2 : Int
}

///|
/// Create a new AlignPath.
pub fn AlignPath::new(
  operations : Array[AlignOp],
  score : Double,
  start1 : Int,
  end1 : Int,
  start2 : Int,
  end2 : Int,
) -> AlignPath {
  AlignPath::{ operations, score, start1, end1, start2, end2 }
}

///|
/// Get the score.
pub fn AlignPath::get_score(self : AlignPath) -> Double {
  self.score
}

///|
/// Get operations.
pub fn AlignPath::get_operations(self : AlignPath) -> Array[AlignOp] {
  self.operations
}

///|
/// Nucleotide substitution matrix (match=+5, mismatch=-4).
fn nucleotide_substitution_matrix(c1 : Char, c2 : Char) -> Double {
  let c1_upper = c1.to_int().to_uint16()
  let c2_upper = c2.to_int().to_uint16()
  let a = 'A'.to_int().to_uint16()
  let t = 'T'.to_int().to_uint16()
  let u = 'U'.to_int().to_uint16()
  let g = 'G'.to_int().to_uint16()
  let c = 'C'.to_int().to_uint16()

  let normalize = fn(c : UInt16) -> UInt16 {
    if c >= 97 && c <= 122 {
      c - 32
    } else {
      c
    }
  }

  let nc1 = normalize(c1_upper)
  let nc2 = normalize(c2_upper)

  if nc1 == nc2 {
    5.0
  } else if (nc1 == a && nc2 == t) ||
    (nc1 == t && nc2 == a) ||
    (nc1 == g && nc2 == c) ||
    (nc1 == c && nc2 == g) ||
    (nc1 == a && nc2 == u) ||
    (nc1 == u && nc2 == a) {
    -4.0
  } else {
    -4.0
  }
}

///|
/// Amino acid characters for BLOSUM62 indexing.
fn amino_acids() -> FixedArray[Char] {
  let arr : FixedArray[Char] = FixedArray::make(20, ' ')
  arr[0] = 'A'
  arr[1] = 'R'
  arr[2] = 'N'
  arr[3] = 'D'
  arr[4] = 'C'
  arr[5] = 'Q'
  arr[6] = 'E'
  arr[7] = 'G'
  arr[8] = 'H'
  arr[9] = 'I'
  arr[10] = 'L'
  arr[11] = 'K'
  arr[12] = 'M'
  arr[13] = 'F'
  arr[14] = 'P'
  arr[15] = 'S'
  arr[16] = 'T'
  arr[17] = 'W'
  arr[18] = 'Y'
  arr[19] = 'V'
  arr
}

///|
/// Get amino acid index (0-19) or -1 for unknown.
fn aa_to_idx(c : Char) -> Int {
  let aa = amino_acids()
  for i = 0; i < 20; i = i + 1 {
    let a = aa[i]
    if c == a || c == a.to_string().to_lower().unsafe_get(0).unsafe_to_char() {
      return i
    }
  }
  -1
}

///|
/// BLOSUM62 substitution matrix lookup.
fn blosum62(c1 : Char, c2 : Char) -> Double {
  let i = aa_to_idx(c1)
  let j = aa_to_idx(c2)
  if i < 0 || j < 0 {
    return -4.0
  }

  let blosum = blosum62_matrix()
  blosum[i][j]
}

///|
/// BLOSUM62 substitution matrix.
fn blosum62_matrix() -> FixedArray[FixedArray[Double]] {
  let mat : FixedArray[FixedArray[Double]] = FixedArray::make(
    20,
    FixedArray::make(20, 0.0),
  )

  mat[0] = make_blosum_row([
    4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0,
    -1.0, 1.0, 0.0, -3.0, -2.0, 0.0,
  ])
  mat[1] = make_blosum_row([
    -1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0,
    -2.0, -1.0, -1.0, -3.0, -2.0, -3.0,
  ])
  mat[2] = make_blosum_row([
    -2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0,
    1.0, 0.0, -4.0, -2.0, -3.0,
  ])
  mat[3] = make_blosum_row([
    -2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0,
    -1.0, 0.0, -1.0, -4.0, -3.0, -3.0,
  ])
  mat[4] = make_blosum_row([
    0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0,
    -3.0, -1.0, -1.0, -2.0, -2.0, -1.0,
  ])
  mat[5] = make_blosum_row([
    -1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0,
    0.0, -1.0, -2.0, -1.0, -2.0,
  ])
  mat[6] = make_blosum_row([
    -1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0,
    -1.0, 0.0, -1.0, -3.0, -2.0, -3.0,
  ])
  mat[7] = make_blosum_row([
    0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0,
    -2.0, 0.0, -2.0, -2.0, -3.0, -3.0,
  ])
  mat[8] = make_blosum_row([
    -2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0,
    -2.0, -1.0, -2.0, -2.0, 2.0, -3.0,
  ])
  mat[9] = make_blosum_row([
    -1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0,
    -3.0, -2.0, -1.0, -3.0, -1.0, 3.0,
  ])
  mat[10] = make_blosum_row([
    -1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0,
    -3.0, -2.0, -1.0, -2.0, -1.0, 1.0,
  ])
  mat[11] = make_blosum_row([
    -1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0,
    -1.0, 0.0, -1.0, -3.0, -2.0, -3.0,
  ])
  mat[12] = make_blosum_row([
    -1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0,
    -2.0, -1.0, -1.0, -1.0, -1.0, 1.0,
  ])
  mat[13] = make_blosum_row([
    -2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0,
    -4.0, -2.0, -2.0, 1.0, 3.0, -1.0,
  ])
  mat[14] = make_blosum_row([
    -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0,
    -4.0, 7.0, -1.0, -1.0, -4.0, -3.0, -2.0,
  ])
  mat[15] = make_blosum_row([
    1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0,
    -1.0, 4.0, 1.0, -3.0, -2.0, -2.0,
  ])
  mat[16] = make_blosum_row([
    0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0,
    -1.0, 1.0, 4.0, -2.0, -2.0, 0.0,
  ])
  mat[17] = make_blosum_row([
    -3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0,
    1.0, -4.0, -3.0, -2.0, 11.0, 2.0, -3.0,
  ])
  mat[18] = make_blosum_row([
    -2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0,
    3.0, -3.0, -2.0, -2.0, 2.0, 7.0, -1.0,
  ])
  mat[19] = make_blosum_row([
    0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -3.0, -3.0, -3.0, 3.0, 1.0, -3.0, 1.0, -1.0,
    -2.0, -2.0, 0.0, -3.0, -1.0, 4.0,
  ])

  mat
}

///|
/// Helper function to create a BLOSUM row.
fn make_blosum_row(values : Array[Double]) -> FixedArray[Double] {
  let row : FixedArray[Double] = FixedArray::make(20, 0.0)
  for i = 0; i < 20; i = i + 1 {
    row[i] = values[i]
  }
  row
}

///|
/// Create string from a slice of FixedArray[Char].
fn fixedarray_char_slice_to_string(
  arr : FixedArray[Char],
  start : Int,
  end : Int,
) -> String {
  let buf = StringBuilder::new(size_hint=end - start)
  for i = start; i < end; i = i + 1 {
    buf.write_char(arr[i])
  }
  buf.to_string()
}

///|
/// Needleman-Wunsch global alignment algorithm.
pub fn global_pairwise_align_nucleotide(
  seq1 : DNA,
  seq2 : DNA,
  gap_open_penalty? : Double = 10.0,
  gap_extend_penalty? : Double = 1.0,
  match_score? : Double = 5.0,
  mismatch_score? : Double = -4.0,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) raise AlignError {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()

  let n = s1.length()
  let m = s2.length()

  let chars1 : FixedArray[Char] = FixedArray::make(n, ' ')
  let chars2 : FixedArray[Char] = FixedArray::make(m, ' ')
  for i = 0; i < n; i = i + 1 {
    chars1[i] = s1.unsafe_get(i).unsafe_to_char()
  }
  for i = 0; i < m; i = i + 1 {
    chars2[i] = s2.unsafe_get(i).unsafe_to_char()
  }

  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0.0),
  )
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0),
  )

  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }

  for i = 1; i <= n; i = i + 1 {
    score_matrix[i][0] = -(gap_open_penalty +
      (i - 1).to_double() * gap_extend_penalty)
    trace_matrix[i][0] = 1
  }

  for j = 1; j <= m; j = j + 1 {
    score_matrix[0][j] = -(gap_open_penalty +
      (j - 1).to_double() * gap_extend_penalty)
    trace_matrix[0][j] = 2
  }

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = chars1[i - 1]
      let c2 = chars2[j - 1]

      let match_score_val = if c1 == c2 { match_score } else { mismatch_score }
      let diag_score = score_matrix[i - 1][j - 1] + match_score_val

      let gap1_score = if trace_matrix[i - 1][j] == 1 {
        score_matrix[i - 1][j] - gap_extend_penalty
      } else {
        score_matrix[i - 1][j] - gap_open_penalty
      }

      let gap2_score = if trace_matrix[i][j - 1] == 2 {
        score_matrix[i][j - 1] - gap_extend_penalty
      } else {
        score_matrix[i][j - 1] - gap_open_penalty
      }

      let max_score = diag_score.max(gap1_score).max(gap2_score)
      score_matrix[i][j] = max_score

      if max_score == diag_score {
        trace_matrix[i][j] = 0
      } else if max_score == gap1_score {
        trace_matrix[i][j] = 1
      } else {
        trace_matrix[i][j] = 2
      }
    }
  }

  let max_len = n + m
  let aligned1_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let aligned2_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let mut pos = max_len - 1

  let mut i = n
  let mut j = m

  while i > 0 || j > 0 {
    match trace_matrix[i][j] {
      0 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = chars2[j - 1]
        i = i - 1
        j = j - 1
      }
      1 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = '-'
        i = i - 1
      }
      2 => {
        aligned1_chars[pos] = '-'
        aligned2_chars[pos] = chars2[j - 1]
        j = j - 1
      }
      _ => break
    }
    pos = pos - 1
  }

  let aligned_seq1 = Seq::new(
    fixedarray_char_slice_to_string(aligned1_chars, pos + 1, max_len),
  )
  let aligned_seq2 = Seq::new(
    fixedarray_char_slice_to_string(aligned2_chars, pos + 1, max_len),
  )

  let record1 = SeqRecord::new(aligned_seq1, id="seq1")
  let record2 = SeqRecord::new(aligned_seq2, id="seq2")

  let records = Array::new()
  records.push(record1)
  records.push(record2)
  let msa = TabularMSA::new(records, dtype="DNA")
  (msa, score_matrix[n][m], (0, n), (0, m))
}

///|
/// Smith-Waterman local alignment algorithm.
pub fn local_pairwise_align_nucleotide(
  seq1 : DNA,
  seq2 : DNA,
  gap_open_penalty? : Double = 10.0,
  gap_extend_penalty? : Double = 1.0,
  match_score? : Double = 5.0,
  mismatch_score? : Double = -4.0,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) raise AlignError {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()

  let n = s1.length()
  let m = s2.length()

  let chars1 : FixedArray[Char] = FixedArray::make(n, ' ')
  let chars2 : FixedArray[Char] = FixedArray::make(m, ' ')
  for i = 0; i < n; i = i + 1 {
    chars1[i] = s1.unsafe_get(i).unsafe_to_char()
  }
  for i = 0; i < m; i = i + 1 {
    chars2[i] = s2.unsafe_get(i).unsafe_to_char()
  }

  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0.0),
  )
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0),
  )

  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }

  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = chars1[i - 1]
      let c2 = chars2[j - 1]

      let match_score_val = if c1 == c2 { match_score } else { mismatch_score }
      let diag_score = score_matrix[i - 1][j - 1] + match_score_val

      let gap1_score = if trace_matrix[i - 1][j] == 1 {
        score_matrix[i - 1][j] - gap_extend_penalty
      } else {
        score_matrix[i - 1][j] - gap_open_penalty
      }

      let gap2_score = if trace_matrix[i][j - 1] == 2 {
        score_matrix[i][j - 1] - gap_extend_penalty
      } else {
        score_matrix[i][j - 1] - gap_open_penalty
      }

      let current_score = diag_score.max(gap1_score).max(gap2_score).max(0.0)
      score_matrix[i][j] = current_score

      if current_score > max_score {
        max_score = current_score
        max_i = i
        max_j = j
      }

      if current_score == 0.0 {
        trace_matrix[i][j] = 0
      } else if current_score == diag_score {
        trace_matrix[i][j] = 1
      } else if current_score == gap1_score {
        trace_matrix[i][j] = 2
      } else {
        trace_matrix[i][j] = 3
      }
    }
  }

  let max_len = n + m
  let aligned1_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let aligned2_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let mut pos = max_len - 1

  let mut i = max_i
  let mut j = max_j

  while i > 0 && j > 0 && score_matrix[i][j] > 0.0 {
    match trace_matrix[i][j] {
      1 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = chars2[j - 1]
        i = i - 1
        j = j - 1
      }
      2 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = '-'
        i = i - 1
      }
      3 => {
        aligned1_chars[pos] = '-'
        aligned2_chars[pos] = chars2[j - 1]
        j = j - 1
      }
      _ => break
    }
    pos = pos - 1
  }

  let aligned_seq1 = Seq::new(
    fixedarray_char_slice_to_string(aligned1_chars, pos + 1, max_len),
  )
  let aligned_seq2 = Seq::new(
    fixedarray_char_slice_to_string(aligned2_chars, pos + 1, max_len),
  )

  let record1 = SeqRecord::new(aligned_seq1, id="seq1")
  let record2 = SeqRecord::new(aligned_seq2, id="seq2")

  let records = Array::new()
  records.push(record1)
  records.push(record2)
  let msa = TabularMSA::new(records, dtype="DNA")
  (msa, max_score, (i, max_i), (j, max_j))
}

///|
/// Global pairwise alignment for protein sequences.
pub fn global_pairwise_align_protein(
  seq1 : Protein,
  seq2 : Protein,
  gap_open_penalty? : Double = 10.0,
  gap_extend_penalty? : Double = 1.0,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) raise AlignError {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()

  let n = s1.length()
  let m = s2.length()

  let chars1 : FixedArray[Char] = FixedArray::make(n, ' ')
  let chars2 : FixedArray[Char] = FixedArray::make(m, ' ')
  for i = 0; i < n; i = i + 1 {
    chars1[i] = s1.unsafe_get(i).unsafe_to_char()
  }
  for i = 0; i < m; i = i + 1 {
    chars2[i] = s2.unsafe_get(i).unsafe_to_char()
  }

  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0.0),
  )
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0),
  )

  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }

  for i = 1; i <= n; i = i + 1 {
    score_matrix[i][0] = -(gap_open_penalty +
      (i - 1).to_double() * gap_extend_penalty)
    trace_matrix[i][0] = 1
  }

  for j = 1; j <= m; j = j + 1 {
    score_matrix[0][j] = -(gap_open_penalty +
      (j - 1).to_double() * gap_extend_penalty)
    trace_matrix[0][j] = 2
  }

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = chars1[i - 1]
      let c2 = chars2[j - 1]

      let diag_score = score_matrix[i - 1][j - 1] + blosum62(c1, c2)

      let gap1_score = if trace_matrix[i - 1][j] == 1 {
        score_matrix[i - 1][j] - gap_extend_penalty
      } else {
        score_matrix[i - 1][j] - gap_open_penalty
      }

      let gap2_score = if trace_matrix[i][j - 1] == 2 {
        score_matrix[i][j - 1] - gap_extend_penalty
      } else {
        score_matrix[i][j - 1] - gap_open_penalty
      }

      let max_score = diag_score.max(gap1_score).max(gap2_score)
      score_matrix[i][j] = max_score

      if max_score == diag_score {
        trace_matrix[i][j] = 0
      } else if max_score == gap1_score {
        trace_matrix[i][j] = 1
      } else {
        trace_matrix[i][j] = 2
      }
    }
  }

  let max_len = n + m
  let aligned1_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let aligned2_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let mut pos = max_len - 1

  let mut i = n
  let mut j = m

  while i > 0 || j > 0 {
    match trace_matrix[i][j] {
      0 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = chars2[j - 1]
        i = i - 1
        j = j - 1
      }
      1 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = '-'
        i = i - 1
      }
      2 => {
        aligned1_chars[pos] = '-'
        aligned2_chars[pos] = chars2[j - 1]
        j = j - 1
      }
      _ => break
    }
    pos = pos - 1
  }

  let aligned_seq1 = Seq::new(
    fixedarray_char_slice_to_string(aligned1_chars, pos + 1, max_len),
  )
  let aligned_seq2 = Seq::new(
    fixedarray_char_slice_to_string(aligned2_chars, pos + 1, max_len),
  )

  let record1 = SeqRecord::new(aligned_seq1, id="seq1")
  let record2 = SeqRecord::new(aligned_seq2, id="seq2")

  let records = Array::new()
  records.push(record1)
  records.push(record2)
  let msa = TabularMSA::new(records, dtype="Protein")
  (msa, score_matrix[n][m], (0, n), (0, m))
}

///|
/// Local pairwise alignment for protein sequences.
pub fn local_pairwise_align_protein(
  seq1 : Protein,
  seq2 : Protein,
  gap_open_penalty? : Double = 10.0,
  gap_extend_penalty? : Double = 1.0,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) raise AlignError {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()

  let n = s1.length()
  let m = s2.length()

  let chars1 : FixedArray[Char] = FixedArray::make(n, ' ')
  let chars2 : FixedArray[Char] = FixedArray::make(m, ' ')
  for i = 0; i < n; i = i + 1 {
    chars1[i] = s1.unsafe_get(i).unsafe_to_char()
  }
  for i = 0; i < m; i = i + 1 {
    chars2[i] = s2.unsafe_get(i).unsafe_to_char()
  }

  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0.0),
  )
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(
    n + 1,
    FixedArray::make(m + 1, 0),
  )

  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }

  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = chars1[i - 1]
      let c2 = chars2[j - 1]

      let diag_score = score_matrix[i - 1][j - 1] + blosum62(c1, c2)

      let gap1_score = if trace_matrix[i - 1][j] == 1 {
        score_matrix[i - 1][j] - gap_extend_penalty
      } else {
        score_matrix[i - 1][j] - gap_open_penalty
      }

      let gap2_score = if trace_matrix[i][j - 1] == 2 {
        score_matrix[i][j - 1] - gap_extend_penalty
      } else {
        score_matrix[i][j - 1] - gap_open_penalty
      }

      let current_score = diag_score.max(gap1_score).max(gap2_score).max(0.0)
      score_matrix[i][j] = current_score

      if current_score > max_score {
        max_score = current_score
        max_i = i
        max_j = j
      }

      if current_score == 0.0 {
        trace_matrix[i][j] = 0
      } else if current_score == diag_score {
        trace_matrix[i][j] = 1
      } else if current_score == gap1_score {
        trace_matrix[i][j] = 2
      } else {
        trace_matrix[i][j] = 3
      }
    }
  }

  let max_len = n + m
  let aligned1_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let aligned2_chars : FixedArray[Char] = FixedArray::make(max_len, ' ')
  let mut pos = max_len - 1

  let mut i = max_i
  let mut j = max_j

  while i > 0 && j > 0 && score_matrix[i][j] > 0.0 {
    match trace_matrix[i][j] {
      1 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = chars2[j - 1]
        i = i - 1
        j = j - 1
      }
      2 => {
        aligned1_chars[pos] = chars1[i - 1]
        aligned2_chars[pos] = '-'
        i = i - 1
      }
      3 => {
        aligned1_chars[pos] = '-'
        aligned2_chars[pos] = chars2[j - 1]
        j = j - 1
      }
      _ => break
    }
    pos = pos - 1
  }

  let aligned_seq1 = Seq::new(
    fixedarray_char_slice_to_string(aligned1_chars, pos + 1, max_len),
  )
  let aligned_seq2 = Seq::new(
    fixedarray_char_slice_to_string(aligned2_chars, pos + 1, max_len),
  )

  let record1 = SeqRecord::new(aligned_seq1, id="seq1")
  let record2 = SeqRecord::new(aligned_seq2, id="seq2")

  let records = Array::new()
  records.push(record1)
  records.push(record2)
  let msa = TabularMSA::new(records, dtype="Protein")
  (msa, max_score, (i, max_i), (j, max_j))
}

///|
/// Compute alignment score.
pub fn align_score(msa : TabularMSA) -> Double {
  let n_rows = msa.num_records()
  let n_cols = msa.get_alignment_length()

  if n_rows < 2 || n_cols == 0 {
    return 0.0
  }

  let mut total_score = 0.0
  let mut count = 0

  for col = 0; col < n_cols; col = col + 1 {
    for i = 0; i < n_rows; i = i + 1 {
      for j = i + 1; j < n_rows; j = j + 1 {
        let c1 = msa.get(i).unwrap().seq.get_char(col).unwrap_or('-')
        let c2 = msa.get(j).unwrap().seq.get_char(col).unwrap_or('-')

        let score = if msa.dtype == "Protein" {
          blosum62(c1, c2)
        } else if c1 == c2 {
          5.0
        } else {
          -4.0
        }

        total_score = total_score + score
        count = count + 1
      }
    }
  }

  if count > 0 {
    total_score / count.to_double()
  } else {
    0.0
  }
}
