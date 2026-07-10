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
    if c >= 97 && c <= 122 { c - 32 }
    else { c }
  }
  
  let nc1 = normalize(c1_upper)
  let nc2 = normalize(c2_upper)
  
  if nc1 == nc2 {
    5.0
  } else if (nc1 == a && nc2 == t) || (nc1 == t && nc2 == a) ||
            (nc1 == g && nc2 == c) || (nc1 == c && nc2 == g) ||
            (nc1 == a && nc2 == u) || (nc1 == u && nc2 == a) {
    -4.0
  } else {
    -4.0
  }
}

///|
/// BLOSUM62 substitution matrix.
fn blosum62(c1 : Char, c2 : Char) -> Double {
  let blosum = Map([
    ("A", Map([("A", 4.0), ("R", -1.0), ("N", -2.0), ("D", -2.0), ("C", 0.0), ("Q", -1.0), ("E", -1.0), ("G", 0.0), ("H", -2.0), ("I", -1.0), ("L", -1.0), ("K", -1.0), ("M", -1.0), ("F", -2.0), ("P", -1.0), ("S", 1.0), ("T", 0.0), ("W", -3.0), ("Y", -2.0), ("V", 0.0)])),
    ("R", Map([("A", -1.0), ("R", 5.0), ("N", 0.0), ("D", -2.0), ("C", -3.0), ("Q", 1.0), ("E", 0.0), ("G", -2.0), ("H", 0.0), ("I", -3.0), ("L", -2.0), ("K", 2.0), ("M", -1.0), ("F", -3.0), ("P", -2.0), ("S", -1.0), ("T", -1.0), ("W", -3.0), ("Y", -2.0), ("V", -3.0)])),
    ("N", Map([("A", -2.0), ("R", 0.0), ("N", 6.0), ("D", 1.0), ("C", -3.0), ("Q", 0.0), ("E", 0.0), ("G", 0.0), ("H", 1.0), ("I", -3.0), ("L", -3.0), ("K", 0.0), ("M", -2.0), ("F", -3.0), ("P", -2.0), ("S", 1.0), ("T", 0.0), ("W", -4.0), ("Y", -2.0), ("V", -3.0)])),
    ("D", Map([("A", -2.0), ("R", -2.0), ("N", 1.0), ("D", 6.0), ("C", -3.0), ("Q", 0.0), ("E", 2.0), ("G", -1.0), ("H", -1.0), ("I", -3.0), ("L", -4.0), ("K", -1.0), ("M", -3.0), ("F", -3.0), ("P", -1.0), ("S", 0.0), ("T", -1.0), ("W", -4.0), ("Y", -3.0), ("V", -3.0)])),
    ("C", Map([("A", 0.0), ("R", -3.0), ("N", -3.0), ("D", -3.0), ("C", 9.0), ("Q", -3.0), ("E", -4.0), ("G", -3.0), ("H", -3.0), ("I", -1.0), ("L", -1.0), ("K", -3.0), ("M", -1.0), ("F", -2.0), ("P", -3.0), ("S", -1.0), ("T", -1.0), ("W", -2.0), ("Y", -2.0), ("V", -1.0)])),
    ("Q", Map([("A", -1.0), ("R", 1.0), ("N", 0.0), ("D", 0.0), ("C", -3.0), ("Q", 5.0), ("E", 2.0), ("G", -2.0), ("H", 0.0), ("I", -3.0), ("L", -2.0), ("K", 1.0), ("M", 0.0), ("F", -3.0), ("P", -1.0), ("S", 0.0), ("T", -1.0), ("W", -2.0), ("Y", -1.0), ("V", -2.0)])),
    ("E", Map([("A", -1.0), ("R", 0.0), ("N", 0.0), ("D", 2.0), ("C", -4.0), ("Q", 2.0), ("E", 5.0), ("G", -2.0), ("H", 0.0), ("I", -3.0), ("L", -3.0), ("K", 1.0), ("M", -2.0), ("F", -3.0), ("P", -1.0), ("S", 0.0), ("T", -1.0), ("W", -3.0), ("Y", -2.0), ("V", -3.0)])),
    ("G", Map([("A", 0.0), ("R", -2.0), ("N", 0.0), ("D", -1.0), ("C", -3.0), ("Q", -2.0), ("E", -2.0), ("G", 6.0), ("H", -2.0), ("I", -4.0), ("L", -4.0), ("K", -2.0), ("M", -3.0), ("F", -3.0), ("P", -2.0), ("S", 0.0), ("T", -2.0), ("W", -2.0), ("Y", -3.0), ("V", -3.0)])),
    ("H", Map([("A", -2.0), ("R", 0.0), ("N", 1.0), ("D", -1.0), ("C", -3.0), ("Q", 0.0), ("E", 0.0), ("G", -2.0), ("H", 8.0), ("I", -3.0), ("L", -3.0), ("K", -1.0), ("M", -2.0), ("F", -1.0), ("P", -2.0), ("S", -1.0), ("T", -2.0), ("W", -2.0), ("Y", 2.0), ("V", -3.0)])),
    ("I", Map([("A", -1.0), ("R", -3.0), ("N", -3.0), ("D", -3.0), ("C", -1.0), ("Q", -3.0), ("E", -3.0), ("G", -4.0), ("H", -3.0), ("I", 4.0), ("L", 2.0), ("K", -3.0), ("M", 1.0), ("F", 0.0), ("P", -3.0), ("S", -2.0), ("T", -1.0), ("W", -3.0), ("Y", -1.0), ("V", 3.0)])),
    ("L", Map([("A", -1.0), ("R", -2.0), ("N", -3.0), ("D", -4.0), ("C", -1.0), ("Q", -2.0), ("E", -3.0), ("G", -4.0), ("H", -3.0), ("I", 2.0), ("L", 4.0), ("K", -2.0), ("M", 2.0), ("F", 0.0), ("P", -3.0), ("S", -2.0), ("T", -1.0), ("W", -2.0), ("Y", -1.0), ("V", 1.0)])),
    ("K", Map([("A", -1.0), ("R", 2.0), ("N", 0.0), ("D", -1.0), ("C", -3.0), ("Q", 1.0), ("E", 1.0), ("G", -2.0), ("H", -1.0), ("I", -3.0), ("L", -2.0), ("K", 5.0), ("M", -1.0), ("F", -3.0), ("P", -1.0), ("S", 0.0), ("T", -1.0), ("W", -3.0), ("Y", -2.0), ("V", -3.0)])),
    ("M", Map([("A", -1.0), ("R", -1.0), ("N", -2.0), ("D", -3.0), ("C", -1.0), ("Q", 0.0), ("E", -2.0), ("G", -3.0), ("H", -2.0), ("I", 1.0), ("L", 2.0), ("K", -1.0), ("M", 5.0), ("F", 0.0), ("P", -2.0), ("S", -1.0), ("T", -1.0), ("W", -1.0), ("Y", -1.0), ("V", 1.0)])),
    ("F", Map([("A", -2.0), ("R", -3.0), ("N", -3.0), ("D", -3.0), ("C", -2.0), ("Q", -3.0), ("E", -3.0), ("G", -3.0), ("H", -1.0), ("I", 0.0), ("L", 0.0), ("K", -3.0), ("M", 0.0), ("F", 6.0), ("P", -4.0), ("S", -2.0), ("T", -2.0), ("W", 1.0), ("Y", 3.0), ("V", -1.0)])),
    ("P", Map([("A", -1.0), ("R", -2.0), ("N", -2.0), ("D", -1.0), ("C", -3.0), ("Q", -1.0), ("E", -1.0), ("G", -2.0), ("H", -2.0), ("I", -3.0), ("L", -3.0), ("K", -1.0), ("M", -2.0), ("F", -4.0), ("P", 7.0), ("S", -1.0), ("T", -1.0), ("W", -4.0), ("Y", -3.0), ("V", -2.0)])),
    ("S", Map([("A", 1.0), ("R", -1.0), ("N", 1.0), ("D", 0.0), ("C", -1.0), ("Q", 0.0), ("E", 0.0), ("G", 0.0), ("H", -1.0), ("I", -2.0), ("L", -2.0), ("K", 0.0), ("M", -1.0), ("F", -2.0), ("P", -1.0), ("S", 4.0), ("T", 1.0), ("W", -3.0), ("Y", -2.0), ("V", -2.0)])),
    ("T", Map([("A", 0.0), ("R", -1.0), ("N", 0.0), ("D", -1.0), ("C", -1.0), ("Q", -1.0), ("E", -1.0), ("G", -2.0), ("H", -2.0), ("I", -1.0), ("L", -1.0), ("K", -1.0), ("M", -1.0), ("F", -2.0), ("P", -1.0), ("S", 1.0), ("T", 4.0), ("W", -2.0), ("Y", -2.0), ("V", 0.0)])),
    ("W", Map([("A", -3.0), ("R", -3.0), ("N", -4.0), ("D", -4.0), ("C", -2.0), ("Q", -2.0), ("E", -3.0), ("G", -2.0), ("H", -2.0), ("I", -3.0), ("L", -2.0), ("K", -3.0), ("M", -1.0), ("F", 1.0), ("P", -4.0), ("S", -3.0), ("T", -2.0), ("W", 11.0), ("Y", 2.0), ("V", -3.0)])),
    ("Y", Map([("A", -2.0), ("R", -2.0), ("N", -2.0), ("D", -3.0), ("C", -2.0), ("Q", -1.0), ("E", -2.0), ("G", -3.0), ("H", 2.0), ("I", -1.0), ("L", -1.0), ("K", -2.0), ("M", -1.0), ("F", 3.0), ("P", -3.0), ("S", -2.0), ("T", -2.0), ("W", 2.0), ("Y", 7.0), ("V", -1.0)])),
    ("V", Map([("A", 0.0), ("R", -3.0), ("N", -3.0), ("D", -3.0), ("C", -1.0), ("Q", -2.0), ("E", -3.0), ("G", -3.0), ("H", -3.0), ("I", 3.0), ("L", 1.0), ("K", -3.0), ("M", 1.0), ("F", -1.0), ("P", -2.0), ("S", -2.0), ("T", 0.0), ("W", -3.0), ("Y", -1.0), ("V", 4.0)])),
  ], capacity=20)
  
  let c1_str = c1.to_string().to_upper()
  let c2_str = c2.to_string().to_upper()
  
  match blosum.get(c1_str) {
    Some(row) => match row.get(c2_str) {
      Some(v) => v
      None => -4.0
    }
    None => -4.0
  }
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
  
  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0.0))
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0))
  
  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }
  
  for i = 1; i <= n; i = i + 1 {
    score_matrix[i][0] = -(gap_open_penalty + (i - 1).to_double() * gap_extend_penalty)
    trace_matrix[i][0] = 1
  }
  
  for j = 1; j <= m; j = j + 1 {
    score_matrix[0][j] = -(gap_open_penalty + (j - 1).to_double() * gap_extend_penalty)
    trace_matrix[0][j] = 2
  }
  
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      
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
  
  let mut i = n
  let mut j = m
  let aligned1 = StringBuilder::new()
  let aligned2 = StringBuilder::new()
  
  while i > 0 || j > 0 {
    match trace_matrix[i][j] {
      0 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      1 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char('-')
        i = i - 1
      }
      2 => {
        aligned1.write_char('-')
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      _ => break
    }
  }
  
  let aligned_seq1 = Seq::new(aligned1.to_string()).reverse()
  let aligned_seq2 = Seq::new(aligned2.to_string()).reverse()
  
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
  
  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0.0))
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0))
  
  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }
  
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0
  
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      
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
  
  let mut i = max_i
  let mut j = max_j
  let aligned1 = StringBuilder::new()
  let aligned2 = StringBuilder::new()
  
  while i > 0 && j > 0 && score_matrix[i][j] > 0.0 {
    match trace_matrix[i][j] {
      1 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      2 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char('-')
        i = i - 1
      }
      3 => {
        aligned1.write_char('-')
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      _ => break
    }
  }
  
  let aligned_seq1 = Seq::new(aligned1.to_string()).reverse()
  let aligned_seq2 = Seq::new(aligned2.to_string()).reverse()
  
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
  
  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0.0))
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0))
  
  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }
  
  for i = 1; i <= n; i = i + 1 {
    score_matrix[i][0] = -(gap_open_penalty + (i - 1).to_double() * gap_extend_penalty)
    trace_matrix[i][0] = 1
  }
  
  for j = 1; j <= m; j = j + 1 {
    score_matrix[0][j] = -(gap_open_penalty + (j - 1).to_double() * gap_extend_penalty)
    trace_matrix[0][j] = 2
  }
  
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      
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
  
  let mut i = n
  let mut j = m
  let aligned1 = StringBuilder::new()
  let aligned2 = StringBuilder::new()
  
  while i > 0 || j > 0 {
    match trace_matrix[i][j] {
      0 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      1 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char('-')
        i = i - 1
      }
      2 => {
        aligned1.write_char('-')
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      _ => break
    }
  }
  
  let aligned_seq1 = Seq::new(aligned1.to_string()).reverse()
  let aligned_seq2 = Seq::new(aligned2.to_string()).reverse()
  
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
  
  let score_matrix : FixedArray[FixedArray[Double]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0.0))
  let trace_matrix : FixedArray[FixedArray[Int]] = FixedArray::make(n + 1, FixedArray::make(m + 1, 0))
  
  for i = 0; i <= n; i = i + 1 {
    score_matrix[i] = FixedArray::make(m + 1, 0.0)
    trace_matrix[i] = FixedArray::make(m + 1, 0)
  }
  
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0
  
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      
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
  
  let mut i = max_i
  let mut j = max_j
  let aligned1 = StringBuilder::new()
  let aligned2 = StringBuilder::new()
  
  while i > 0 && j > 0 && score_matrix[i][j] > 0.0 {
    match trace_matrix[i][j] {
      1 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      2 => {
        aligned1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
        aligned2.write_char('-')
        i = i - 1
      }
      3 => {
        aligned1.write_char('-')
        aligned2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      _ => break
    }
  }
  
  let aligned_seq1 = Seq::new(aligned1.to_string()).reverse()
  let aligned_seq2 = Seq::new(aligned2.to_string()).reverse()
  
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
        } else {
          if c1 == c2 { 5.0 } else { -4.0 }
        }
        
        total_score = total_score + score
        count = count + 1
      }
    }
  }
  
  if count > 0 { total_score / count.to_double() } else { 0.0 }
}