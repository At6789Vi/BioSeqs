///|
/// CRAM format parsing module.
/// CRAM is a compressed binary format for sequence alignment data.

///|
/// CRAM parsing errors.
pub suberror CramError {
  CramError(String)
}

///|
/// CRAM magic number ("CRAM").
let cram_magic : Array[Int] = [67, 82, 65, 77]

///|
/// Check if data starts with CRAM magic number.
pub fn is_cram_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  for i = 0; i < 4; i = i + 1 {
    if data.unsafe_get(i).to_int() != cram_magic[i] {
      return false
    }
  }
  true
}

///|
/// CRAM file definition structure.
pub struct CramFileDefinition {
  major_version : Int
  minor_version : Int
  flags : Int
}

///|
/// CRAM block compression methods.
pub enum CramCompressionMethod {
  NoCompression
  Gzip
  Bzip2
  Lzma
  Rans
  Bgzf
}

///|
/// CRAM block content types.
pub enum CramContentType {
  CompressionHeader
  SliceHeader
  CoreData
  ExternalData
  FileHeader
  ReferenceSequence
  Mapping
  Info
  QualityScore
  TagData
  ReadFeatures
  InsertSize
  SubstitutionMatrix
  DeletionMatrix
  InsertionMatrix
  Padding
  Unknown(Int)
}

///|
/// CRAM block structure.
pub struct CramBlock {
  content_type : CramContentType
  compression_method : CramCompressionMethod
  compressed_size : Int
  raw_size : Int
  data : String
}

///|
/// CRAM slice header.
pub struct CramSliceHeader {
  alignment_start : Int
  alignment_span : Int
  num_records : Int
  num_blocks : Int
  reference_id : Int
  landmark : Int
}

///|
/// CRAM compression header.
pub struct CramCompressionHeader {
  version : Int
  reference_id : Int
  num_ref_bases : Int
  ref_seq_id : Int
  ref_start : Int
  ref_end : Int
  num_slices : Int
}

///|
/// CRAM read feature types.
pub enum CramReadFeature {
  Substitution(Int, Int)
  Insertion(Int, String)
  Deletion(Int, Int)
  SkippedBases(Int, Int)
  HardClip(Int, Int)
  SoftClip(Int, String)
  Padding(Int, Int)
  QualityScore(Int, String)
  ModifiedBase(Int, String)
  ModifiedQuality(Int, String)
  UnknownFeature(Int)
}

///|
/// CRAM alignment record.
pub struct CramRecord {
  seq_id : Int
  alignment_start : Int
  read_length : Int
  read_group : Int
  flags : Int
  mapping_quality : Int
  mate_seq_id : Int
  mate_alignment_start : Int
  template_size : Int
  num_reads : Int
  read_number : Int
  fragment_length : Int
  qname : String
  cigar : FixedArray[BamCigarElement]
  seq : String
  qual : String
  tags : Array[(String, String)]
  read_features : FixedArray[CramReadFeature]
}

///|
/// CRAM file representation.
pub struct CramFile {
  header : BamHeader
  records : Array[CramRecord]
}

///|
/// Convert content type code to enum.
fn content_type_from_code(code : Int) -> CramContentType {
  match code {
    0 => CramContentType::CompressionHeader
    1 => CramContentType::SliceHeader
    2 => CramContentType::CoreData
    3 => CramContentType::ExternalData
    4 => CramContentType::FileHeader
    5 => CramContentType::ReferenceSequence
    6 => CramContentType::Mapping
    7 => CramContentType::Info
    8 => CramContentType::QualityScore
    9 => CramContentType::TagData
    10 => CramContentType::ReadFeatures
    11 => CramContentType::InsertSize
    12 => CramContentType::SubstitutionMatrix
    13 => CramContentType::DeletionMatrix
    14 => CramContentType::InsertionMatrix
    255 => CramContentType::Padding
    _ => CramContentType::Unknown(code)
  }
}

///|
/// Convert compression method code to enum.
fn compression_method_from_code(code : Int) -> CramCompressionMethod {
  match code {
    0 => CramCompressionMethod::NoCompression
    1 => CramCompressionMethod::Gzip
    2 => CramCompressionMethod::Bzip2
    3 => CramCompressionMethod::Lzma
    4 => CramCompressionMethod::Rans
    5 => CramCompressionMethod::Bgzf
    _ => CramCompressionMethod::NoCompression
  }
}

///|
/// Read ITF-8 integer (Integer Transformation Format - 8 bit).
fn read_itf8(data : String, pos : Int) -> (Int, Int) {
  let mut result = 0
  let mut shift = 0
  let mut current_pos = pos
  let byte = data.unsafe_get(current_pos).to_int() & 0xFF

  let prefix = byte >> 3
  let bytes_to_read = match prefix {
    0 => 1
    1 => 2
    3 => 3
    7 => 4
    15 => 5
    _ => 1
  }

  current_pos = current_pos + 1

  if bytes_to_read == 1 {
    result = byte & 0x7F
  } else {
    result = byte & 0x7
    let mut j = 1
    while j < bytes_to_read {
      let b = data.unsafe_get(current_pos).to_int() & 0xFF
      result = result | (b << shift)
      shift = shift + 8
      current_pos = current_pos + 1
      j = j + 1
    }
  }

  if (result & 1) == 1 {
    result = -1 - (result >> 1)
  } else {
    result = result >> 1
  }

  (result, current_pos)
}

///|
/// Read LTF-8 long integer.
fn read_ltf8(data : String, pos : Int) -> (Int, Int) {
  let byte = data.unsafe_get(pos).to_int() & 0xFF
  let prefix = byte >> 3
  let bytes_to_read = match prefix {
    0 => 1
    1 => 2
    3 => 3
    7 => 4
    15 => 5
    31 => 6
    63 => 7
    127 => 8
    _ => 1
  }

  let mut result = 0
  let mut shift = 0
  let mut current_pos = pos + 1

  if bytes_to_read == 1 {
    result = byte & 0x7F
  } else {
    result = byte & 0x7
    let mut j = 1
    while j < bytes_to_read {
      let b = data.unsafe_get(current_pos).to_int() & 0xFF
      result = result | (b << shift)
      shift = shift + 8
      current_pos = current_pos + 1
      j = j + 1
    }
  }

  if (result & 1) == 1 {
    result = -1 - (result >> 1)
  } else {
    result = result >> 1
  }

  (result, current_pos)
}

///|
/// Parse CRAM file definition.
fn parse_cram_file_definition(
  data : String,
  pos : Int,
) -> (CramFileDefinition, Int) {
  let (major_version, p1) = read_int8(data, pos + 4)
  let (minor_version, p2) = read_int8(data, p1)
  let (flags, p3) = read_int8(data, p2)

  (CramFileDefinition::{ major_version, minor_version, flags }, p3 + 5)
}

///|
/// Parse CRAM block.
fn parse_cram_block(data : String, pos : Int) -> (CramBlock, Int) {
  let (content_type_code, p1) = read_int8(data, pos)
  let (compression_code, p2) = read_int8(data, p1)
  let (compressed_size, p3) = read_uint32_le(data, p2)
  let (raw_size, p4) = read_uint32_le(data, p3)

  let content = data[p4:p4 + compressed_size].to_owned()

  (
    CramBlock::{
      content_type: content_type_from_code(content_type_code),
      compression_method: compression_method_from_code(compression_code),
      compressed_size,
      raw_size,
      data: content,
    },
    p4 + compressed_size,
  )
}

///|
/// Parse CRAM compression header block.
fn parse_compression_header(data : String) -> CramCompressionHeader {
  let mut pos = 0
  let (version, p1) = read_itf8(data, pos)
  pos = p1

  let (reference_id, p2) = read_itf8(data, pos)
  pos = p2

  let (num_ref_bases, p3) = read_itf8(data, pos)
  pos = p3

  let (ref_seq_id, p4) = read_itf8(data, pos)
  pos = p4

  let (ref_start, p5) = read_itf8(data, pos)
  pos = p5

  let (ref_end, p6) = read_itf8(data, pos)
  pos = p6

  let (num_slices, p7) = read_itf8(data, pos)
  pos = p7

  CramCompressionHeader::{
    version,
    reference_id,
    num_ref_bases,
    ref_seq_id,
    ref_start,
    ref_end,
    num_slices,
  }
}

///|
/// Parse CRAM slice header block.
fn parse_slice_header(data : String) -> CramSliceHeader {
  let mut pos = 0
  let (alignment_start, p1) = read_itf8(data, pos)
  pos = p1

  let (alignment_span, p2) = read_itf8(data, pos)
  pos = p2

  let (num_records, p3) = read_itf8(data, pos)
  pos = p3

  let (num_blocks, p4) = read_itf8(data, pos)
  pos = p4

  let (reference_id, p5) = read_itf8(data, pos)
  pos = p5

  let (landmark, p6) = read_itf8(data, pos)
  pos = p6

  CramSliceHeader::{
    alignment_start,
    alignment_span,
    num_records,
    num_blocks,
    reference_id,
    landmark,
  }
}

///|
/// Parse CRAM read features.
fn parse_read_features(
  data : String,
  length : Int,
) -> FixedArray[CramReadFeature] {
  let features : FixedArray[CramReadFeature] = FixedArray::make(
    length,
    CramReadFeature::UnknownFeature(0),
  )
  let mut pos = 0
  let mut idx = 0

  while pos < data.length() && idx < length {
    let (code, p1) = read_int8(data, pos)
    pos = p1

    match code {
      1 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (base, p3) = read_int8(data, pos)
        pos = p3
        features[idx] = CramReadFeature::Substitution(pos_offset, base)
      }
      2 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (insert_len, p3) = read_itf8(data, pos)
        pos = p3
        let insert_seq = data[pos:pos + insert_len].to_owned()
        pos = pos + insert_len
        features[idx] = CramReadFeature::Insertion(pos_offset, insert_seq)
      }
      3 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (del_len, p3) = read_itf8(data, pos)
        pos = p3
        features[idx] = CramReadFeature::Deletion(pos_offset, del_len)
      }
      4 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (skip_len, p3) = read_itf8(data, pos)
        pos = p3
        features[idx] = CramReadFeature::SkippedBases(pos_offset, skip_len)
      }
      5 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (clip_len, p3) = read_itf8(data, pos)
        pos = p3
        features[idx] = CramReadFeature::HardClip(pos_offset, clip_len)
      }
      6 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (clip_len, p3) = read_itf8(data, pos)
        pos = p3
        let clip_seq = data[pos:pos + clip_len].to_owned()
        pos = pos + clip_len
        features[idx] = CramReadFeature::SoftClip(pos_offset, clip_seq)
      }
      7 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (pad_len, p3) = read_itf8(data, pos)
        pos = p3
        features[idx] = CramReadFeature::Padding(pos_offset, pad_len)
      }
      8 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (qual_len, p3) = read_itf8(data, pos)
        pos = p3
        let qual_str = data[pos:pos + qual_len].to_owned()
        pos = pos + qual_len
        features[idx] = CramReadFeature::QualityScore(pos_offset, qual_str)
      }
      9 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (mod_len, p3) = read_itf8(data, pos)
        pos = p3
        let mod_str = data[pos:pos + mod_len].to_owned()
        pos = pos + mod_len
        features[idx] = CramReadFeature::ModifiedBase(pos_offset, mod_str)
      }
      10 => {
        let (pos_offset, p2) = read_itf8(data, pos)
        pos = p2
        let (mod_len, p3) = read_itf8(data, pos)
        pos = p3
        let mod_str = data[pos:pos + mod_len].to_owned()
        pos = pos + mod_len
        features[idx] = CramReadFeature::ModifiedQuality(pos_offset, mod_str)
      }
      _ => features[idx] = CramReadFeature::UnknownFeature(code)
    }

    idx = idx + 1
  }

  features
}

///|
/// Parse CRAM core data block and extract records.
fn parse_cram_core_data(data : String, header : BamHeader) -> Array[CramRecord] {
  let records = Array::new()
  let mut pos = 0

  while pos < data.length() {
    let (seq_id, p1) = read_itf8(data, pos)
    pos = p1

    let (alignment_start, p2) = read_itf8(data, pos)
    pos = p2

    let (read_length, p3) = read_itf8(data, pos)
    pos = p3

    let (read_group, p4) = read_itf8(data, pos)
    pos = p4

    let (flags, p5) = read_itf8(data, pos)
    pos = p5

    let (mapping_quality, p6) = read_itf8(data, pos)
    pos = p6

    let (mate_seq_id, p7) = read_itf8(data, pos)
    pos = p7

    let (mate_alignment_start, p8) = read_itf8(data, pos)
    pos = p8

    let (template_size, p9) = read_itf8(data, pos)
    pos = p9

    let (num_reads, p10) = read_itf8(data, pos)
    pos = p10

    let (read_number, p11) = read_itf8(data, pos)
    pos = p11

    let (fragment_length, p12) = read_itf8(data, pos)
    pos = p12

    let (qname, p13) = read_cstring(data, pos)
    pos = p13

    let (cigar_len, p14) = read_itf8(data, pos)
    pos = p14

    let cigar : FixedArray[BamCigarElement] = FixedArray::make(cigar_len, BamCigarElement::{
      op: BamCigarOp::Match,
      length: 0,
    })
    let mut cigar_idx = 0
    while cigar_idx < cigar_len && pos < data.length() {
      let (cigar_op, p15) = read_itf8(data, pos)
      pos = p15
      let (cigar_len_val, p16) = read_itf8(data, pos)
      pos = p16
      cigar[cigar_idx] = BamCigarElement::{
        op: cigar_op_from_code(cigar_op),
        length: cigar_len_val,
      }
      cigar_idx = cigar_idx + 1
    }

    let (seq_len, p17) = read_itf8(data, pos)
    pos = p17

    let seq = data[pos:pos + seq_len].to_owned()
    pos = pos + seq_len

    let (qual_len, p18) = read_itf8(data, pos)
    pos = p18

    let qual = data[pos:pos + qual_len].to_owned()
    pos = pos + qual_len

    let (num_tags, p19) = read_itf8(data, pos)
    pos = p19

    let tags = Array::new()
    let mut tag_idx = 0
    while tag_idx < num_tags && pos < data.length() {
      let (tag_name, p20) = read_cstring(data, pos)
      pos = p20
      let (tag_value, p21) = read_cstring(data, pos)
      pos = p21
      tags.push((tag_name, tag_value))
      tag_idx = tag_idx + 1
    }

    let (num_features, p22) = read_itf8(data, pos)
    pos = p22

    let read_features = parse_read_features(data, num_features)
    pos = pos + num_features * 4

    records.push(CramRecord::{
      seq_id,
      alignment_start,
      read_length,
      read_group,
      flags,
      mapping_quality,
      mate_seq_id,
      mate_alignment_start,
      template_size,
      num_reads,
      read_number,
      fragment_length,
      qname,
      cigar,
      seq,
      qual,
      tags,
      read_features,
    })
  }

  records
}

///|
/// Parse raw CRAM data (uncompressed).
pub fn parse_cram(data : String) -> CramFile raise CramError {
  let n = data.length()
  if n < 12 {
    raise CramError("CRAM file too short")
  }

  for i = 0; i < 4; i = i + 1 {
    if data.unsafe_get(i).to_int() != cram_magic[i] {
      raise CramError("Invalid CRAM magic number")
    }
  }

  let (file_def, pos) = parse_cram_file_definition(data, 0)

  let mut current_pos = pos
  let mut bam_header : BamHeader = BamHeader::{
    text: "",
    references: FixedArray::make(0, BamReference::{ name: "", length: 0 }),
  }
  let records = Array::new()

  while current_pos < n {
    let (container_len, p1) = read_int32_le(data, current_pos)
    current_pos = p1

    let (ref_id, p2) = read_itf8(data, current_pos)
    current_pos = p2

    let (blocks_len, p3) = read_uint32_le(data, current_pos)
    current_pos = p3

    let (num_records, p4) = read_itf8(data, current_pos)
    current_pos = p4

    let (record_counter, p5) = read_ltf8(data, current_pos)
    current_pos = p5

    let (base_count, p6) = read_ltf8(data, current_pos)
    current_pos = p6

    let (container_checksum, p7) = read_uint32_le(data, current_pos)
    current_pos = p7

    let (block_count, p8) = read_itf8(data, current_pos)
    current_pos = p8

    let mut block_idx = 0
    while block_idx < block_count && current_pos < n {
      let (block, p9) = parse_cram_block(data, current_pos)
      current_pos = p9

      match block.content_type {
        CramContentType::FileHeader =>
          bam_header = BamHeader::{
            text: block.data,
            references: FixedArray::make(0, BamReference::{
              name: "",
              length: 0,
            }),
          }
        CramContentType::CoreData => {
          let core_records = parse_cram_core_data(block.data, bam_header)
          for rec in core_records {
            records.push(rec)
          }
        }
        _ => ()
      }

      block_idx = block_idx + 1
    }

    current_pos = current_pos + container_len - (current_pos - (pos + 4))
  }

  CramFile::{ header: bam_header, records }
}

///|
/// Get the number of references in the CRAM file.
pub fn CramFile::num_references(self : CramFile) -> Int {
  self.header.references.length()
}

///|
/// Get reference names from the CRAM header.
pub fn CramFile::get_reference_names(self : CramFile) -> Array[String] {
  let names = Array::new()
  for ref_info in self.header.references {
    names.push(ref_info.name)
  }
  names
}

///|
/// Get a reference name by index.
pub fn CramRecord::get_reference_name(
  self : CramRecord,
  header : BamHeader,
) -> String {
  if self.seq_id >= 0 && self.seq_id < header.references.length() {
    header.references[self.seq_id].name
  } else {
    "*"
  }
}

///|
/// Convert CRAM record to SAM-like string representation.
pub fn CramRecord::to_sam_string(
  self : CramRecord,
  header : BamHeader,
) -> String {
  let rname = self.get_reference_name(header)
  let cigar_str = StringBuilder::new()
  for elem in self.cigar {
    let op_char = match elem.op {
      BamCigarOp::Match => 'M'
      BamCigarOp::Insert => 'I'
      BamCigarOp::Delete => 'D'
      BamCigarOp::Skip => 'N'
      BamCigarOp::SoftClip => 'S'
      BamCigarOp::HardClip => 'H'
      BamCigarOp::Pad => 'P'
      BamCigarOp::Equal => '='
      BamCigarOp::Diff => 'X'
    }
    cigar_str.write_string(elem.length.to_string())
    cigar_str.write_char(op_char)
  }

  let pos_str = if self.alignment_start >= 0 {
    (self.alignment_start + 1).to_string()
  } else {
    "0"
  }

  "\{self.qname}\t\{self.flags}\t\{rname}\t\{pos_str}\t\{self.mapping_quality}\t\{cigar_str.to_string()}\t*\t0\t0\t\{self.seq}\t\{self.qual}"
}

///|
/// Convert CRAM file to SAM string.
pub fn CramFile::to_sam(self : CramFile) -> String {
  let buf = StringBuilder::new()

  buf.write_string(self.header.text)

  for record in self.records {
    buf.write_string(record.to_sam_string(self.header))
    buf.write_char('\n')
  }

  buf.to_string()
}

///|
/// Get the end position of the alignment.
pub fn CramRecord::get_end(self : CramRecord) -> Int {
  if self.alignment_start == 0 || self.cigar.length() == 0 {
    return 0
  }

  let mut end = self.alignment_start - 1
  let mut i = 0
  while i < self.cigar.length() {
    let element = self.cigar[i]
    let op = element.op
    if op == BamCigarOp::Match ||
      op == BamCigarOp::Delete ||
      op == BamCigarOp::Skip ||
      op == BamCigarOp::Equal ||
      op == BamCigarOp::Diff {
      end = end + element.length
    }
    i = i + 1
  }
  end
}

///|
/// Check if the read is paired.
pub fn CramRecord::is_paired(self : CramRecord) -> Bool {
  (self.flags & 1) != 0
}

///|
/// Check if the read is mapped.
pub fn CramRecord::is_mapped(self : CramRecord) -> Bool {
  (self.flags & 4) == 0
}

///|
/// Check if the read is reverse complemented.
pub fn CramRecord::is_reverse(self : CramRecord) -> Bool {
  (self.flags & 16) != 0
}
