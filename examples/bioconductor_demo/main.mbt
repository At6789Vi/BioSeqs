///|
fn main {
  println("=== ChIPseeker Peak Annotation Example ===")
  println("Inspired by Bioconductor ChIPseeker package")
  println("")

  let peaks = [
    ("chr1", 1000, 2000),
    ("chr1", 5000, 6000),
    ("chr2", 10000, 11000),
  ]

  let genes = [
    @bio.Gene::new("gene1", "chr1", 1, 5000, @bio.strand_plus()),
    @bio.Gene::new("gene2", "chr2", 10001, 20000, @bio.strand_minus()),
  ]
  let transcripts : Array[@bio.Transcript] = []
  let exons : Array[@bio.Exon] = []
  let cds : Array[@bio.CDS] = []
  let txdb = @bio.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @bio.annotate_peaks(peaks, txdb)

  println("1. Peak Annotation Results")
  println("   Total annotations: \{result.annotations.length().to_string()}")
  println("")

  let mut i = 0
  while i < result.annotations.length() && i < 3 {
    let ann = result.annotations[i]
    println("   Peak \{i + 1}: \{ann.chr}:\{ann.start}-\{ann.end}")
    println("     - Annotation: \{ann.annotation}")
    println("     - Gene: \{ann.gene_name}")
    println("     - Distance to TSS: \{ann.distance_to_tss.to_string()}")
    println("     - Feature type: \{ann.feature_type}")
    println("")
    i = i + 1
  }

  println("2. Annotation Summary")
  let keys = result.summary.keys().collect()
  let mut j = 0
  while j < keys.length() {
    let key = keys[j]
    println("   \{key}: \{result.summary[key].to_string()}")
    j = j + 1
  }
  println("")

  println("3. Annotation Statistics")
  let stats = @bio.get_annotation_statistics(result)
  let stat_keys = stats.keys().collect()
  let mut k = 0
  while k < stat_keys.length() {
    let stat_key = stat_keys[k]
    println("   \{stat_key}: \{stats[stat_key].to_string()}")
    k = k + 1
  }
  println("")

  println("4. Promoter Annotations")
  let promoters = @bio.get_promoter_annotations(result)
  println("   Promoter peaks: \{promoters.length().to_string()}")
  println("")

  println("5. Annotation Plot")
  println(@bio.plot_annotation(result.summary))
  println("")

  println("\n")
  println("=== DOSE Disease Ontology Enrichment Example ===")
  println("Inspired by Bioconductor DOSE package")
  println("")

  let gene_list = ["TP53", "BRCA1", "BRCA2", "ATM", "CDKN1A"]
  let background_genes = [
    "TP53", "BRCA1", "BRCA2", "ATM", "CDKN1A", "GAPDH", "ACTB", "MYC", "EGFR", "VEGFA",
  ]

  let disease_db = Map([], capacity=4)
  disease_db["DOID:1234"] = ["TP53", "BRCA1", "CDKN1A"]
  disease_db["DOID:5678"] = ["BRCA2", "ATM"]
  disease_db["DOID:9012"] = ["GAPDH", "ACTB"]
  disease_db["DOID:3456"] = ["MYC", "EGFR"]

  let dose_result = @bio.enrich_disease(gene_list, disease_db, background_genes)

  println("1. Enrichment Results")
  println("   Total enriched terms: \{dose_result.terms.length().to_string()}")
  println("   Query genes: \{dose_result.query_count.to_string()}")
  println("   Background genes: \{dose_result.background_count.to_string()}")
  println("")

  println("2. Top 3 Enriched Terms")
  let top_terms = @bio.get_top_terms(dose_result, 3)
  i = 0
  while i < top_terms.length() {
    let term = top_terms[i]
    println("   Term \{i + 1}:")
    println("     - ID: \{term.term_id}")
    println("     - Name: \{term.term_name}")
    println("     - Gene count: \{term.gene_count.to_string()}")
    println("     - P-value: \{term.p_value.to_string()}")
    println("     - Adjusted P-value: \{term.adjusted_p_value.to_string()}")
    println("     - Enrichment score: \{term.enrichment_score.to_string()}")
    println("")
    i = i + 1
  }

  println("3. Filtered Results (adjusted p-value < 0.05)")
  let filtered = @bio.filter_results(dose_result, 0.05)
  println("   Significant terms: \{filtered.terms.length().to_string()}")
  println("")

  println("4. Summary")
  println(@bio.summarize_result(dose_result))
  println("")

  println("\n")
  println("=== ReactomePA Pathway Analysis Example ===")
  println("Inspired by Bioconductor ReactomePA package")
  println("")

  let reactome_gene_list = ["TP53", "BRCA1", "ATM", "CDKN1A"]
  let reactome_bg_genes = [
    "TP53", "BRCA1", "ATM", "CDKN1A", "GAPDH", "ACTB", "MYC",
  ]

  let reactome_db = Map([], capacity=4)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1", "CDKN1A"],
    "Cellular Processes",
  )
  reactome_db["R-HSA-5678"] = (
    "DNA Repair",
    ["BRCA1", "ATM"],
    "DNA Replication",
  )
  reactome_db["R-HSA-9012"] = (
    "Metabolism",
    ["GAPDH", "ACTB"],
    "Metabolic Processes",
  )
  reactome_db["R-HSA-3456"] = (
    "Signal Transduction",
    ["MYC"],
    "Signal Transduction",
  )

  let reactome_result = @bio.enrich_reactome(
    reactome_gene_list, reactome_db, reactome_bg_genes,
  )

  println("1. Pathway Enrichment Results")
  println(
    "   Total enriched pathways: \{reactome_result.pathways.length().to_string()}",
  )
  println("   Query genes: \{reactome_result.query_count.to_string()}")
  println(
    "   Background genes: \{reactome_result.background_count.to_string()}",
  )
  println("")

  println("2. Top 3 Pathways")
  let top_pathways = @bio.get_top_pathways(reactome_result, 3)
  i = 0
  while i < top_pathways.length() {
    let pathway = top_pathways[i]
    println("   Pathway \{i + 1}:")
    println("     - ID: \{pathway.pathway_id}")
    println("     - Name: \{pathway.pathway_name}")
    println("     - Gene count: \{pathway.gene_count.to_string()}")
    println("     - P-value: \{pathway.p_value.to_string()}")
    println("     - Adjusted P-value: \{pathway.adjusted_p_value.to_string()}")
    println("     - Enrichment score: \{pathway.enrichment_score.to_string()}")
    println("     - Top level term: \{pathway.top_level_term}")
    println("")
    i = i + 1
  }

  println("3. Pathways Grouped by Top Level Term")
  let groups = @bio.group_by_top_level(reactome_result)
  let group_keys = groups.keys().collect()
  j = 0
  while j < group_keys.length() {
    let group_key = group_keys[j]
    println(
      "   \{group_key}: \{groups[group_key].length().to_string()} pathways",
    )
    j = j + 1
  }
  println("")

  println("4. Filtered Results (adjusted p-value < 0.05)")
  let reactome_filtered = @bio.reactome_filter_results(reactome_result, 0.05)
  println(
    "   Significant pathways: \{reactome_filtered.pathways.length().to_string()}",
  )
  println("")

  println("5. Summary")
  println(@bio.reactome_summarize_result(reactome_result))
}
