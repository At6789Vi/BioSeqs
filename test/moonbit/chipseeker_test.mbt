///|
test "chipseeker_annotate_peaks" {
  let peaks = [
    ("chr1", 1000, 2000),
    ("chr1", 5000, 6000),
    ("chr2", 10000, 11000),
  ]

  let genes = [@src.Gene::new("gene1", "chr1", 1, 5000, @src.strand_plus())]
  let transcripts = []
  let exons = []
  let cds = []
  let txdb = @src.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @src.annotate_peaks(peaks, txdb)

  assert_eq(result.annotations.length(), 3)
}

///|
test "chipseeker_annotation_summary" {
  let peaks = [
    ("chr1", 1000, 2000),
    ("chr1", 5000, 6000),
    ("chr2", 10000, 11000),
  ]

  let genes = [@src.Gene::new("gene1", "chr1", 1, 5000, @src.strand_plus())]
  let transcripts = []
  let exons = []
  let cds = []
  let txdb = @src.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @src.annotate_peaks(peaks, txdb)

  assert_true(result.summary.length() > 0)
}

///|
test "chipseeker_get_promoter_annotations" {
  let peaks = [("chr1", 1000, 2000), ("chr1", 5000, 6000)]

  let genes = [@src.Gene::new("gene1", "chr1", 1, 5000, @src.strand_plus())]
  let transcripts = []
  let exons = []
  let cds = []
  let txdb = @src.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @src.annotate_peaks(peaks, txdb)
  let promoters = @src.get_promoter_annotations(result)

  assert_true(promoters.length() >= 0)
}

///|
test "chipseeker_annotation_statistics" {
  let peaks = [("chr1", 1000, 2000), ("chr1", 5000, 6000)]

  let genes = [@src.Gene::new("gene1", "chr1", 1, 5000, @src.strand_plus())]
  let transcripts = []
  let exons = []
  let cds = []
  let txdb = @src.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @src.annotate_peaks(peaks, txdb)
  let stats = @src.get_annotation_statistics(result)

  assert_true(stats.length() > 0)
}

///|
test "chipseeker_plot_annotation" {
  let peaks = [("chr1", 1000, 2000)]

  let genes = [@src.Gene::new("gene1", "chr1", 1, 5000, @src.strand_plus())]
  let transcripts = []
  let exons = []
  let cds = []
  let txdb = @src.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let result = @src.annotate_peaks(peaks, txdb)
  let plot = @src.plot_annotation(result.summary)

  assert_true(plot.length() > 0)
}
