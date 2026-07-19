///|
fn main {
  let data = @bio.gsva_create_example_data()
  let gene_sets = @bio.gsva_create_example_gene_sets()

  let params = @bio.GSVAParams::new()
  let scores = @bio.gsva_run(data, gene_sets, params)

  println(
    "GSVA Scores computed for " +
    scores.scores.length().to_string() +
    " gene sets",
  )
  println("Number of samples: " + scores.sample_names.length().to_string())

  let summary = @bio.gsva_summarize(scores)
  println(summary)
}
