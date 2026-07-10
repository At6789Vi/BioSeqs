///|
test "crc32_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = crc32(seq)
  let unsigned = if result >= 0 { result } else { result + 4294967296 }
  assert_eq(unsigned, 4247645206)
}

///|
test "gcg_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = gcg(seq)
  assert_eq(result, 5493)
}

///|
test "seguid_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = seguid(seq)
  assert_eq(result, "kspO3Jbmvk9xX3jVv8QOSapllIo")
}

///|
test "gc_simple" {
  let result = gc("AGCT")
  assert_true(result > 49.9 && result < 50.1)
}

///|
test "gc_full" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = gc(seq)
  assert_true(result > 41.0 && result < 42.0)
}

///|
test "protein_molecular_weight" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.molecular_weight()
  assert_true(result > 6380.0 && result < 6390.0)
}

///|
test "protein_aromaticity" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.aromaticity()
  assert_true(result > 0.1 && result < 0.11)
}

///|
test "protein_gravy" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.gravy()
  assert_true(result > -0.2 && result < -0.1)
}

///|
test "tm_wallace_simple" {
  let result = tm_wallace("AAAAATTTTT")
  assert_eq(result, 20.0)
}

///|
test "tm_wallace_gc" {
  let result = tm_wallace("GGCCGGCC")
  assert_eq(result, 32.0)
}

///|
test "tm_gc_simple" {
  let result = tm_gc("AGCTAGCT")
  assert_true(result > 5.0 && result < 6.0)
}

///|
test "tm_gc_long" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = tm_gc(seq)
  assert_true(result > 68.0 && result < 69.0)
}