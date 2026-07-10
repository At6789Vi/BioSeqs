///|
/// Tests for scikit-bio alignment module.

///|
test "dna_creation" {
  let dna = DNA::new("ACGT")
  assert_eq(dna.to_string(), "ACGT")
  assert_eq(dna.length(), 4)
}

///|
test "rna_creation" {
  let rna = RNA::new("ACGU")
  assert_eq(rna.to_string(), "ACGU")
  assert_eq(rna.length(), 4)
}

///|
test "protein_creation" {
  let prot = Protein::new("ACDEFGHIKLMNPQRSTVWY")
  assert_eq(prot.to_string(), "ACDEFGHIKLMNPQRSTVWY")
  assert_eq(prot.length(), 20)
}

///|
test "tabular_msa_creation" {
  let dna1 = DNA::new("ACGT")
  let dna2 = DNA::new("ACCT")
  let seqs : Array[DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = TabularMSA::from_dna(seqs)

  assert_eq(msa.num_records(), 2)
  assert_eq(msa.get_alignment_length(), 4)
}

///|
test "tabular_msa_column" {
  let dna1 = DNA::new("ACGT")
  let dna2 = DNA::new("ACCT")
  let seqs : Array[DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = TabularMSA::from_dna(seqs)

  assert_eq(msa.column(0).unwrap(), "AA")
  assert_eq(msa.column(1).unwrap(), "CC")
  assert_eq(msa.column(2).unwrap(), "GC")
  assert_eq(msa.column(3).unwrap(), "TT")
}

///|
test "tabular_msa_consensus" {
  let dna1 = DNA::new("ACGT")
  let dna2 = DNA::new("ACCT")
  let dna3 = DNA::new("AGGT")
  let seqs : Array[DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  seqs.push(dna3)
  let msa = TabularMSA::from_dna(seqs)

  let consensus = msa.consensus()
  assert_eq(consensus.to_string(), "ACGT")
}

///|
test "tabular_msa_conservation" {
  let dna1 = DNA::new("ACGT")
  let dna2 = DNA::new("ACCT")
  let dna3 = DNA::new("AGGT")
  let seqs : Array[DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  seqs.push(dna3)
  let msa = TabularMSA::from_dna(seqs)

  let col2 = msa.column(2).unwrap()
  assert_eq(col2, "GCG")

  let conservation = msa.conservation()
  assert_eq(conservation.length(), 4)
  assert_eq(conservation[0], 1.0)
  assert_eq(conservation[1], 0.6666666666666666)
  assert_eq(conservation[2], 0.6666666666666666)
  assert_eq(conservation[3], 1.0)
}

///|
test "global_pairwise_align_nucleotide" {
  let seq1 = DNA::new("ACGT")
  let seq2 = DNA::new("AGT")

  let (msa, score, pos1, pos2) = global_pairwise_align_nucleotide(seq1, seq2)

  assert_eq(msa.num_records(), 2)
  assert_eq(pos1.0, 0)
  assert_eq(pos1.1, 4)
  assert_eq(pos2.0, 0)
  assert_eq(pos2.1, 3)
}

///|
test "local_pairwise_align_nucleotide" {
  let seq1 = DNA::new("ACGTACGT")
  let seq2 = DNA::new("CGT")

  let (msa, score, pos1, pos2) = local_pairwise_align_nucleotide(seq1, seq2)

  assert_eq(msa.num_records(), 2)

  let col0 = msa.column(0).unwrap()
  assert_eq(col0, "CC")
  let col1 = msa.column(1).unwrap()
  assert_eq(col1, "GG")
  let col2 = msa.column(2).unwrap()
  assert_eq(col2, "TT")
}

///|
test "global_pairwise_align_protein" {
  let seq1 = Protein::new("ACDE")
  let seq2 = Protein::new("ACE")

  let (msa, score, pos1, pos2) = global_pairwise_align_protein(seq1, seq2)

  assert_eq(msa.num_records(), 2)
}

///|
test "local_pairwise_align_protein" {
  let seq1 = Protein::new("ACDEFGHI")
  let seq2 = Protein::new("EFG")

  let (msa, score, pos1, pos2) = local_pairwise_align_protein(seq1, seq2)

  assert_eq(msa.num_records(), 2)
}

///|
test "align_score" {
  let dna1 = DNA::new("ACGT")
  let dna2 = DNA::new("ACGT")
  let seqs : Array[DNA] = Array::new()
  seqs.push(dna1)
  seqs.push(dna2)
  let msa = TabularMSA::from_dna(seqs)

  let score = align_score(msa)
  assert_eq(score, 5.0)
}

///|
test "global_alignment_score" {
  let seq1 = DNA::new("ATCG")
  let seq2 = DNA::new("ATAG")

  let (msa, _, _, _) = global_pairwise_align_nucleotide(seq1, seq2)

  assert_eq(msa.get_alignment_length(), 4)
}
