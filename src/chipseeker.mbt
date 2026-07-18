///|
pub struct PeakAnnotation {
  peak_id : String
  chr : String
  start : Int
  end : Int
  annotation : String
  gene_name : String
  distance_to_tss : Int
  feature_type : String
}

///|
pub struct AnnotationResult {
  annotations : Array[PeakAnnotation]
  summary : Map[String, Int]
}

///|
pub fn annotate_peaks(
  peaks : Array[(String, Int, Int)],
  txdb : TxDb,
) -> AnnotationResult {
  let annotations : Array[PeakAnnotation] = []
  let summary = Map([], capacity=16)

  let mut i = 0
  while i < peaks.length() {
    let (chr, start, end) = peaks[i]
    let peak_id = "peak_${(i + 1).to_string()}"

    let annotation = annotate_single_peak(chr, start, end, txdb)
    annotations.push(annotation)

    let feature_type = annotation.feature_type
    let count = if summary.contains(feature_type) {
      summary[feature_type]
    } else {
      0
    }
    summary[feature_type] = count + 1

    i = i + 1
  }

  AnnotationResult::{ annotations, summary }
}

///|
fn annotate_single_peak(
  chr : String,
  start : Int,
  end : Int,
  txdb : TxDb,
) -> PeakAnnotation {
  let peak_center = (start + end) / 2

  let genes = gene_ids(txdb)

  if genes.length() == 0 {
    return PeakAnnotation::{
      peak_id: "",
      chr,
      start,
      end,
      annotation: "Intergenic",
      gene_name: "",
      distance_to_tss: 1000000,
      feature_type: "Intergenic",
    }
  }

  let mut min_distance = 1000000
  let mut closest_gene = genes[0]

  let mut j = 0
  while j < genes.length() {
    let gene = genes[j]
    let gene_start = 1
    let gene_end = 10000

    let distance = if peak_center < gene_start {
      gene_start - peak_center
    } else if peak_center > gene_end {
      peak_center - gene_end
    } else {
      0
    }

    if distance < min_distance {
      min_distance = distance
      closest_gene = gene
    }
    j = j + 1
  }

  let (annotation, feature_type) = classify_annotation(min_distance)

  PeakAnnotation::{
    peak_id: "",
    chr,
    start,
    end,
    annotation,
    gene_name: closest_gene,
    distance_to_tss: min_distance,
    feature_type,
  }
}

///|
fn classify_annotation(distance : Int) -> (String, String) {
  if distance == 0 {
    ("Promoter", "Promoter")
  } else if distance <= 1000 {
    ("Promoter (<=1kb)", "Promoter")
  } else if distance <= 3000 {
    ("Promoter (1-3kb)", "Promoter")
  } else if distance <= 5000 {
    ("Promoter (3-5kb)", "Promoter")
  } else if distance <= 10000 {
    ("Distal Intergenic", "Distal")
  } else {
    ("Intergenic", "Intergenic")
  }
}

///|
pub fn plot_annotation(summary : Map[String, Int]) -> String {
  let mut result = "Annotation Distribution:\n"
  let keys = summary.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    result = result + "- "
    result = result + key
    result = result + ": "
    result = result + summary[key].to_string()
    result = result + "\n"
    i = i + 1
  }
  result
}

///|
pub fn get_annotation_statistics(
  result : AnnotationResult,
) -> Map[String, Double] {
  let stats = Map([], capacity=16)
  let total = result.annotations.length().to_double()

  let keys = result.summary.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let key = keys[i]
    let value = result.summary[key]
    let percentage = value.to_double() / total * 100.0
    stats[key] = percentage
    i = i + 1
  }

  stats
}

///|
pub fn get_promoter_annotations(
  result : AnnotationResult,
) -> Array[PeakAnnotation] {
  let promoters : Array[PeakAnnotation] = []
  let mut i = 0
  while i < result.annotations.length() {
    let annotation = result.annotations[i]
    if annotation.feature_type == "Promoter" {
      promoters.push(annotation)
    }
    i = i + 1
  }
  promoters
}
