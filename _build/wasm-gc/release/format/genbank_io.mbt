// GenBank format parser, mirroring BioPython's Bio.SeqIO.GenBankIO.
//
// GenBank is a complex multi-section format:
//
//   LOCUS       NM_001     5028 bp    DNA    linear   PRI 21-JUN-1999
//   DEFINITION  Description text here.
//   ACCESSION   NM_001
//   VERSION     NM_001.1  GI:12345
//   SOURCE      Organism name
//   FEATURES    ...
//   ORIGIN
//           1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//          61 ...
//   //
//
// Key fields extracted:
//   - id:          VERSION (first token, e.g. "NM_001.1")
//   - name:        LOCUS name (first token after LOCUS)
//   - description: DEFINITION text (joined, trailing '.' stripped)
//   - seq:         ORIGIN section (numbers and spaces stripped)

///|
/// Parse GenBank format text and return an array of SeqRecord objects.
///
/// Mirrors BioPython's SeqIO.parse(handle, "genbank").
/// Extracts: id (VERSION), name (LOCUS), description (DEFINITION),
/// and sequence (ORIGIN section). Multiple records separated by '//'.
pub fn parse_genbank(content : String) -> Array[SeqRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i].trim_end(chars="\r").to_owned()
    // Look for LOCUS line to start a record
    if line.has_prefix("LOCUS") {
      let mut name = "<unknown>"
      let mut id = "<unknown id>"
      let mut description = "<unknown description>"
      let seq_parts = StringBuilder::new()
      let mut in_origin = false
      // Parse LOCUS line: "LOCUS   NAME   LENGTH bp ..."
      let locus_rest = line[5:].trim_start(chars=" ").to_owned()
      name = match locus_rest.find(" ") {
        Some(pos) => locus_rest[0:pos].to_owned()
        None => locus_rest
      }
      i = i + 1
      // Parse subsequent lines until //
      while i < n {
        let l = lines[i].trim_end(chars="\r").to_owned()
        if l == "//" {
          break
        }
        if l.has_prefix("DEFINITION") {
          // Join continuation lines (indented)
          let def_part = l[10:].trim_start(chars=" ").to_owned()
          let def_buf = StringBuilder::new()
          def_buf.write_string(def_part)
          i = i + 1
          while i < n {
            let cont = lines[i].trim_end(chars="\r").to_owned()
            if cont.length() == 0 ||
              cont.unsafe_get(0) != ' '.to_int().to_uint16() {
              break
            }
            def_buf.write_string(" ")
            def_buf.write_string(cont.trim(chars=" ").to_owned())
            i = i + 1
          }
          description = def_buf.to_string()
          // BioPython strips trailing dot
          if description.has_prefix(".") {
            // shouldn't start with .
          }
          if description.length() > 0 &&
            description.unsafe_get(description.length() - 1) ==
            '.'.to_int().to_uint16() {
            description = description[0:description.length() - 1].to_owned()
          }
          continue
        }
        if l.has_prefix("VERSION") {
          let ver_part = l[7:].trim_start(chars=" ").to_owned()
          id = match ver_part.find(" ") {
            Some(pos) => ver_part[0:pos].to_owned()
            None => ver_part
          }
        }
        if l.has_prefix("ORIGIN") {
          in_origin = true
          i = i + 1
          continue
        }
        if in_origin {
          // Sequence line: remove numbers and spaces, uppercase letters
          // (BioPython normalizes GenBank sequences to uppercase)
          for c in l.iter() {
            if c >= 'a' && c <= 'z' {
              seq_parts.write_char((c.to_int() - 32).unsafe_to_char())
            } else if c >= 'A' && c <= 'Z' {
              seq_parts.write_char(c)
            }
          }
        }
        i = i + 1
      }
      records.push(SeqRecord::{
        seq: Seq::new(seq_parts.to_string()),
        id,
        name,
        description,
        letter_annotations: Map([], capacity=0),
      })
    }
    i = i + 1
  }
  records
}
