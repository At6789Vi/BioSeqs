///|
/// Suffix Array and Suffix Tree Demo
fn main {
  println("=== Suffix Array and Suffix Tree Demo ===")
  println("")

  println("1. Suffix Array - Basic operations:")
  let sa = @src.SuffixArray::new("banana")
  println("   Text: banana")
  println("   Suffix array length: \{sa.len()}")
  println("   Sorted suffixes:")
  for i = 0; i < sa.len(); i = i + 1 {
    println("     [\{i}] \{sa.get_sorted_suffix(i)}")
  }
  println("")

  println("2. Suffix Array - Pattern matching:")
  println("   Count of 'ana': \{sa.count("ana")}")
  println("   Count of 'an': \{sa.count("an")}")
  println("   Count of 'b': \{sa.count("b")}")
  println("   Count of 'xyz': \{sa.count("xyz")}")
  let pos = sa.locate("ana")
  let pos_str = if pos.length() == 0 {
    "none"
  } else {
    let buf = StringBuilder::new()
    for i = 0; i < pos.length(); i = i + 1 {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(pos[i].to_string())
    }
    buf.to_string()
  }
  println("   Positions of 'ana': \{pos_str}")
  println("")

  println("3. LCP Array:")
  let lcp = @src.LCPArray::new(sa)
  println("   LCP values: \{lcp.lcp()}")
  println("   Max LCP: \{lcp.max_lcp()}")
  println("")

  println("4. Suffix Tree - Basic operations:")
  let st = @src.SuffixTree::new("banana")
  println("   Text: banana")
  println("   Tree size (nodes): \{st.size()}")
  println("   Contains 'ana': \{st.contains("ana")}")
  println("   Contains 'ban': \{st.contains("ban")}")
  println("   Contains 'xyz': \{st.contains("xyz")}")
  println("")

  println("5. Suffix Tree - Pattern matching:")
  println("   Count of 'ana': \{st.count("ana")}")
  println("   Count of 'an': \{st.count("an")}")
  println("   Positions of 'ana': \{st.locate("ana")}")
  println("")

  println("6. Suffix Tree - Longest repeated substring:")
  println("   LRS of 'banana': \{st.longest_repeated_substring()}")
  println("")

  println("7. DNA sequence example:")
  let dna_sa = @src.SuffixArray::new("ACGTACGT")
  let dna_st = @src.SuffixTree::new("ACGTACGT")
  println("   Text: ACGTACGT")
  println("   SA count of 'ACGT': \{dna_sa.count("ACGT")}")
  println("   ST count of 'ACGT': \{dna_st.count("ACGT")}")
  let dna_pos = dna_sa.locate("ACGT")
  let dna_pos_str = if dna_pos.length() == 0 {
    "none"
  } else {
    let buf = StringBuilder::new()
    for i = 0; i < dna_pos.length(); i = i + 1 {
      if i > 0 {
        buf.write_string(", ")
      }
      buf.write_string(dna_pos[i].to_string())
    }
    buf.to_string()
  }
  println("   SA positions of 'ACGT': \{dna_pos_str}")
  println("")

  println("8. Suffix Tree structure:")
  let st_abc = @src.SuffixTree::new("abc")
  println(st_abc.to_string())
}
