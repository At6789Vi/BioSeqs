// FASTQ format parser and writer, mirroring BioPython's
// Bio.SeqIO.QualityIO.
//
// FASTQ format (Sanger, Phred+33):
//   @sequence_id optional_description
//   ACGTACGT...
//   +optional_repeat_of_header
//   !"#$%&'...
//
// Line 1: @ + title (first word = id, whole title = description)
// Line 2: raw sequence
// Line 3: + (optionally followed by the title again)
// Line 4: quality string (same length as sequence)
//
// Quality scores: Q = ord(char) - 33 (Phred+33 / Sanger encoding)

///|
/// Parse FASTQ format text (Sanger/Phred+33) and return an array of
/// SeqRecord objects with phred_quality letter annotations.
///
/// Mirrors BioPython's SeqIO.parse(handle, "fastq").
pub fn parse_fastq(content : String) -> Array[SeqRecord] raise SeqIOError {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  while i + 3 < n {
    let line0 = lines[i].trim_end(chars="\r").to_owned()
    let line1 = lines[i + 1].trim_end(chars="\r").to_owned()
    let line2 = lines[i + 2].trim_end(chars="\r").to_owned()
    let line3 = lines[i + 3].trim_end(chars="\r").to_owned()
    // Skip trailing empty lines
    if line0.length() == 0 && i + 4 >= n {
      break
    }
    if line0.length() == 0 || line0.unsafe_get(0) != '@'.to_int().to_uint16() {
      raise SeqIOError("Expected '@' at start of FASTQ header, got: " + line0)
    }
    if line2.length() == 0 || line2.unsafe_get(0) != '+'.to_int().to_uint16() {
      raise SeqIOError(
        "Expected '+' at start of FASTQ separator, got: " + line2,
      )
    }
    if line1.length() != line3.length() {
      raise SeqIOError(
        "Sequence length " +
        line1.length().to_string() +
        " and quality length " +
        line3.length().to_string() +
        " don't match in FASTQ record",
      )
    }
    // Parse title
    let title = line0[1:].to_owned()
    let (id, description) = match title.find(" ") {
      Some(pos) => (title[0:pos].to_owned(), title)
      None => (title, title)
    }
    // Decode quality scores (Phred+33)
    let quals = Array::new()
    for c in line3.iter() {
      quals.push(c.to_int() - 33)
    }
    let letter_ann : Map[String, Array[Int]] = Map([], capacity=1)
    letter_ann["phred_quality"] = quals
    records.push(SeqRecord::{
      seq: Seq::new(line1),
      id,
      name: id,
      description,
      letter_annotations: letter_ann,
    })
    i = i + 4
  }
  records
}

///|
/// Write an array of SeqRecord objects to FASTQ format text (Sanger/Phred+33).
///
/// Mirrors BioPython's SeqIO.write(records, handle, "fastq").
/// - Records must have phred_quality letter annotations.
/// - Quality scores are encoded as chr(Q + 33).
pub fn write_fastq(records : Array[SeqRecord]) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let quals = match rec.letter_annotations.get("phred_quality") {
      Some(q) => q
      None =>
        raise SeqIOError(
          "Missing phred_quality annotation for record " + rec.id,
        )
    }
    let seq = rec.seq.to_string()
    if seq.length() != quals.length() {
      raise SeqIOError(
        "Sequence length " +
        seq.length().to_string() +
        " and quality length " +
        quals.length().to_string() +
        " don't match in record " +
        rec.id,
      )
    }
    // Title line: "@" + description (matching BioPython)
    buf.write_char('@')
    buf.write_string(build_fasta_title(rec.id, rec.description))
    buf.write_char('\n')
    // Sequence line
    buf.write_string(seq)
    buf.write_char('\n')
    // Separator line
    buf.write_string("+\n")
    // Quality line: encode each Q as chr(Q + 33)
    let qbuf : FixedArray[UInt16] = FixedArray::make(quals.length(), 0)
    for j = 0; j < quals.length(); j = j + 1 {
      qbuf[j] = (quals[j] + 33).to_uint16()
    }
    buf.write_string(unsafe_from_uint16_fixedarray(qbuf))
    buf.write_char('\n')
  }
  buf.to_string()
}
