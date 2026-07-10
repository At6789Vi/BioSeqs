///|
/// Benchmark test for alignment performance.

fn init {
  println("Running global pairwise align nucleotide benchmark (100x100)...")
  let start = System::time()
  for _ = 0; _ < 10; _ = _ + 1 {
    let seq1 = DNA::new(random_dna(100))
    let seq2 = DNA::new(random_dna(100))
    global_pairwise_align_nucleotide(seq1, seq2)
  }
  let end = System::time()
  println("  Time: {} ms", (end - start) * 1000.0)
  
  println("Running global pairwise align nucleotide benchmark (200x200)...")
  let start2 = System::time()
  for _ = 0; _ < 5; _ = _ + 1 {
    let seq1 = DNA::new(random_dna(200))
    let seq2 = DNA::new(random_dna(200))
    global_pairwise_align_nucleotide(seq1, seq2)
  }
  let end2 = System::time()
  println("  Time: {} ms", (end2 - start2) * 1000.0)
  
  println("Running local pairwise align nucleotide benchmark (100x100)...")
  let start3 = System::time()
  for _ = 0; _ < 10; _ = _ + 1 {
    let seq1 = DNA::new(random_dna(100))
    let seq2 = DNA::new(random_dna(100))
    local_pairwise_align_nucleotide(seq1, seq2)
  }
  let end3 = System::time()
  println("  Time: {} ms", (end3 - start3) * 1000.0)
  
  println("Running global pairwise align protein benchmark (50x50)...")
  let start4 = System::time()
  for _ = 0; _ < 10; _ = _ + 1 {
    let seq1 = Protein::new(random_protein(50))
    let seq2 = Protein::new(random_protein(50))
    global_pairwise_align_protein(seq1, seq2)
  }
  let end4 = System::time()
  println("  Time: {} ms", (end4 - start4) * 1000.0)
  
  println("Running local pairwise align protein benchmark (50x50)...")
  let start5 = System::time()
  for _ = 0; _ < 10; _ = _ + 1 {
    let seq1 = Protein::new(random_protein(50))
    let seq2 = Protein::new(random_protein(50))
    local_pairwise_align_protein(seq1, seq2)
  }
  let end5 = System::time()
  println("  Time: {} ms", (end5 - start5) * 1000.0)
  
  println("All benchmarks completed!")
}
