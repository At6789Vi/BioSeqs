// FASTA format parser and writer, mirroring BioPython's
// Bio.SeqIO.FastaIO.
//
// FASTA format:
//   >header_line
//   sequence_line_1
//   sequence_line_2
//   ...
//   >next_header
//   ...
//
// The header line (after '>') is the title. The first whitespace-delimited
// token becomes the record id; the entire title becomes the description.

///|
/// Default line width for wrapping FASTA sequence output (matches BioPython).
let fasta_line_width : Int = 60

///|
/// Parse FASTA format text and return an array of SeqRecord objects.
///
/// Mirrors BioPython's SimpleFastaParser + SeqIO.parse(handle, "fasta").
/// - The title line (after '>') is split at the first whitespace: the first
///   token becomes `id`, the whole title becomes `description`.
/// - `name` is set equal to `id`.
/// - Sequence lines are concatenated (whitespace stripped).
pub fn parse_fasta(content : String) -> Array[SeqRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let mut current_title : String? = None
  let mut current_seq = StringBuilder::new()
  // Walk lines; when we hit a '>' line, flush the previous record.
  for line in lines {
    let view = line.trim_end(chars="\r")
    let trimmed = view.to_owned()
    if trimmed.length() > 0 && trimmed.unsafe_get(0) == '>'.to_int().to_uint16() {
      // Flush previous record
      match current_title {
        Some(title) => {
          let seq_str = current_seq.to_string()
          records.push(fasta_title_to_record(title, seq_str))
        }
        None => ()
      }
      current_title = Some(trimmed[1:].to_owned())
      current_seq = StringBuilder::new()
    } else {
      // Sequence line — strip all whitespace and append
      let seq_part = trimmed.replace(old=" ", new="").replace(old="\t", new="")
      if seq_part.length() > 0 {
        current_seq.write_string(seq_part)
      }
    }
  }
  // Flush last record
  match current_title {
    Some(title) => {
      let seq_str = current_seq.to_string()
      records.push(fasta_title_to_record(title, seq_str))
    }
    None => ()
  }
  records
}

///|
/// Build a SeqRecord from a FASTA title and sequence string.
/// The title is split at the first whitespace: first token = id, whole
/// title = description (matching BioPython).
fn fasta_title_to_record(title : String, seq_str : String) -> SeqRecord {
  let (id, description) = match title.find(" ") {
    Some(pos) => (title[0:pos].to_owned(), title)
    None => (title, title)
  }
  SeqRecord::{
    seq: Seq::new(seq_str),
    id,
    name: id,
    description,
    letter_annotations: Map([], capacity=0),
  }
}

///|
/// Write an array of SeqRecord objects to FASTA format text.
///
/// Mirrors BioPython's FastaIO writer:
/// - The title line is ">" + description (if description starts with id)
///   or ">" + id + " " + description (otherwise).
/// - Sequences are wrapped at 60 characters per line.
pub fn write_fasta(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder::new()
  for rec in records {
    // Determine the title line, matching BioPython's FastaIO logic:
    //   - If description starts with id: title = description
    //   - If description is non-empty but doesn't start with id:
    //     title = id + " " + description
    //   - If description is empty: title = id
    let title = build_fasta_title(rec.id, rec.description)
    buf.write_char('>')
    buf.write_string(title)
    buf.write_char('\n')
    // Write sequence wrapped at fasta_line_width
    let seq = rec.seq.to_string()
    let n = seq.length()
    let mut i = 0
    while i < n {
      let end = if i + fasta_line_width < n { i + fasta_line_width } else { n }
      buf.write_string(seq[i:end].to_owned())
      buf.write_char('\n')
      i = end
    }
  }
  buf.to_string()
}

///|
/// Build the FASTA title line from id and description.
/// If description starts with id, use description as-is.
/// Otherwise, prepend id + " ".
fn build_fasta_title(id : String, description : String) -> String {
  if description.length() == 0 {
    id
  } else if description.has_prefix(id) {
    description
  } else {
    id + " " + description
  }
}
