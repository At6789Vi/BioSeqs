///|
/// PDB format parser and writer for protein structures.

///|
/// Parse a PDB file content into a Structure.
pub fn parse_pdb(content : String) -> Structure {
  let lines = content.split("\n").to_array()
  let models = Array::new()
  let mut current_model_chains = Array::new()
  let mut current_chain_residues = Array::new()
  let mut current_residue_atoms = Array::new()
  let mut current_resname = ""
  let mut current_chainid = ' '
  let mut current_resseq = 1
  let mut current_hetfield = " "
  
  for line in lines {
    let trimmed = line.trim().to_string()
    if trimmed.length() < 6 {
      continue
    }
    
    let record_type = trimmed[0:6].trim().to_string()
    
    if record_type == "MODEL" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(resname=current_resname, chainid=current_chainid, resseq=current_resseq, atoms=current_residue_atoms, hetfield=current_hetfield)
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
      if current_model_chains.length() > 0 {
        let m = Model::new(id=models.length() + 1, chains=current_model_chains)
        models.push(m)
        current_model_chains = Array::new()
      }
    } else if record_type == "ATOM" || record_type == "HETATM" {
      let atom = parse_atom_line(trimmed, record_type == "HETATM")
      
      if current_residue_atoms.length() > 0 && (atom.chainid != current_chainid || atom.resseq != current_resseq) {
        let r = Residue::new(resname=current_resname, chainid=current_chainid, resseq=current_resseq, atoms=current_residue_atoms, hetfield=current_hetfield)
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
        
        if atom.chainid != current_chainid {
          let c = Chain::new(id=current_chainid, residues=current_chain_residues)
          current_model_chains.push(c)
          current_chain_residues = Array::new()
        }
      }
      
      current_resname = atom.resname
      current_chainid = atom.chainid
      current_resseq = atom.resseq
      current_hetfield = if record_type == "HETATM" { "H" } else { " " }
      current_residue_atoms.push(atom)
    } else if record_type == "TER" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(resname=current_resname, chainid=current_chainid, resseq=current_resseq, atoms=current_residue_atoms, hetfield=current_hetfield)
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
    } else if record_type == "ENDMDL" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(resname=current_resname, chainid=current_chainid, resseq=current_resseq, atoms=current_residue_atoms, hetfield=current_hetfield)
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
      if current_model_chains.length() > 0 {
        let m = Model::new(id=models.length() + 1, chains=current_model_chains)
        models.push(m)
        current_model_chains = Array::new()
      }
    }
  }
  
  if current_residue_atoms.length() > 0 {
    let r = Residue::new(resname=current_resname, chainid=current_chainid, resseq=current_resseq, atoms=current_residue_atoms, hetfield=current_hetfield)
    current_chain_residues.push(r)
  }
  if current_chain_residues.length() > 0 {
    let c = Chain::new(id=current_chainid, residues=current_chain_residues)
    current_model_chains.push(c)
  }
  if current_model_chains.length() > 0 {
    let m = Model::new(id=models.length() + 1, chains=current_model_chains)
    models.push(m)
  }
  
  Structure::new(id="", models=models)
}

///|
/// Parse an ATOM or HETATM line into an Atom.
fn parse_atom_line(line : String, _is_hetatm : Bool) -> Atom {
  let name = if line.length() >= 16 { line[12:16].trim().to_string() } else { "" }
  let resname = if line.length() >= 20 { line[17:20].trim().to_string() } else { "UNK" }
  let chainid = if line.length() > 21 { line.get_char(21).unwrap_or(' ') } else { ' ' }
  let resseq_str = if line.length() >= 26 { line[22:26].trim().to_string() } else { "1" }
  let resseq = pdb_parse_int(resseq_str).unwrap_or(1)
  
  let x_str = if line.length() >= 38 { line[30:38].trim().to_string() } else { "0.0" }
  let y_str = if line.length() >= 46 { line[38:46].trim().to_string() } else { "0.0" }
  let z_str = if line.length() >= 54 { line[46:54].trim().to_string() } else { "0.0" }
  let x = pdb_parse_double(x_str).unwrap_or(0.0)
  let y = pdb_parse_double(y_str).unwrap_or(0.0)
  let z = pdb_parse_double(z_str).unwrap_or(0.0)
  
  let occupancy_str = if line.length() >= 60 { line[54:60].trim().to_string() } else { "1.0" }
  let occupancy = pdb_parse_double(occupancy_str).unwrap_or(1.0)
  
  let bfactor_str = if line.length() >= 66 { line[60:66].trim().to_string() } else { "0.0" }
  let bfactor = pdb_parse_double(bfactor_str).unwrap_or(0.0)
  
  let element = if line.length() >= 78 { line[76:78].trim().to_string() } else { "" }
  
  Atom::new(name=name, coord=Vector3::new(x, y, z), resname=resname, chainid=chainid, resseq=resseq, occupancy=occupancy, bfactor=bfactor, element=element, charge=0)
}

///|
/// Parse an integer from a string.
fn pdb_parse_int(s : String) -> Int? {
  let n = s.length()
  if n == 0 {
    None
  } else {
    let mut result = 0
    let zero = '0'.to_int()
    let mut i = 0
    let mut is_negative = false
    
    if s.unsafe_get(0).to_int() == 45 {
      is_negative = true
      i = 1
    }
    
    while i < n {
      let c = s.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 {
        return None
      }
      result = result * 10 + (c - zero)
      i = i + 1
    }
    
    if is_negative {
      Some(-result)
    } else {
      Some(result)
    }
  }
}

///|
/// Parse a double from a string.
fn pdb_parse_double(s : String) -> Double? {
  let n = s.length()
  if n == 0 {
    None
  } else {
    let mut result = 0.0
    let mut decimal = 0.0
    let mut decimal_places = 0
    let mut is_negative = false
    let zero = '0'.to_int()
    let mut i = 0
    
    if s.unsafe_get(0).to_int() == 45 {
      is_negative = true
      i = 1
    }
    
    while i < n {
      let c = s.unsafe_get(i).to_int()
      if c == 46 {
        i = i + 1
        while i < n {
          let dc = s.unsafe_get(i).to_int()
          if dc < zero || dc > zero + 9 {
            break
          }
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
          i = i + 1
        }
        break
      }
      if c < zero || c > zero + 9 {
        return None
      }
      result = result * 10.0 + (c - zero).to_double()
      i = i + 1
    }
    
    while decimal_places > 0 {
      decimal = decimal / 10.0
      decimal_places = decimal_places - 1
    }
    
    let final_result = result + decimal
    if is_negative {
      Some(-final_result)
    } else {
      Some(final_result)
    }
  }
}

///|
/// Write a Structure to PDB format string.
pub fn write_pdb(structure : Structure) -> String {
  let buf = StringBuilder::new()
  
  let mut model_idx = 0
  for model in structure.models {
    model_idx = model_idx + 1
    buf.write_string("MODEL     ")
    buf.write_string(model_idx.to_string())
    buf.write_char('\n')
    
    for chain in model.chains {
      for residue in chain.residues {
        for atom in residue.atoms {
          buf.write_string(write_atom_line(atom, residue.hetfield == "H"))
          buf.write_char('\n')
        }
        buf.write_string("TER\n")
      }
    }
    
    buf.write_string("ENDMDL\n")
  }
  
  buf.write_string("END\n")
  buf.to_string()
}

///|
/// Write an atom to a PDB ATOM or HETATM line.
fn write_atom_line(atom : Atom, is_hetatm : Bool) -> String {
  let buf = StringBuilder::new()
  
  buf.write_string(if is_hetatm { "HETATM" } else { "ATOM  " })
  
  let serial_str = pad_left("1", 5, ' ')
  buf.write_string(serial_str)
  buf.write_char(' ')
  
  let name_str = pad_right(atom.name, 4, ' ')
  buf.write_string(name_str)
  
  buf.write_char(' ')
  buf.write_char(' ')
  
  let resname_str = pad_right(atom.resname, 3, ' ')
  buf.write_string(resname_str)
  buf.write_char(' ')
  
  buf.write_char(atom.chainid)
  buf.write_char(' ')
  
  let resseq_str = pad_left(atom.resseq.to_string(), 4, ' ')
  buf.write_string(resseq_str)
  
  buf.write_char(' ')
  buf.write_string("   ")
  
  let x_str = format_double(atom.coord.x, 8, 3)
  buf.write_string(x_str)
  
  let y_str = format_double(atom.coord.y, 8, 3)
  buf.write_string(y_str)
  
  let z_str = format_double(atom.coord.z, 8, 3)
  buf.write_string(z_str)
  
  let occ_str = format_double(atom.occupancy, 6, 2)
  buf.write_string(occ_str)
  
  let bfactor_str = format_double(atom.bfactor, 6, 2)
  buf.write_string(bfactor_str)
  
  buf.write_string("      ")
  
  let elem_str = pad_right(atom.element, 2, ' ')
  buf.write_string(elem_str)
  
  buf.write_string("  ")
  
  buf.to_string()
}

///|
/// Pad string to the left.
fn pad_left(s : String, width : Int, pad_char : Char) -> String {
  let n = s.length()
  if n >= width {
    s
  } else {
    pad_char.to_string().repeat(width - n) + s
  }
}

///|
/// Pad string to the right.
fn pad_right(s : String, width : Int, pad_char : Char) -> String {
  let n = s.length()
  if n >= width {
    s
  } else {
    s + pad_char.to_string().repeat(width - n)
  }
}

///|
/// Format a double with specified width and precision.
fn format_double(value : Double, width : Int, precision : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_string() } else { "0" }
  let frac_part = if parts.length() > 1 { parts[1].to_string() } else { "" }
  
  let mut padded_frac = frac_part
  if padded_frac.length() < precision {
    padded_frac = padded_frac + "0".repeat(precision - padded_frac.length())
  } else if padded_frac.length() > precision {
    padded_frac = padded_frac[0:precision].to_string()
  }
  
  let result = int_part + "." + padded_frac
  pad_left(result, width, ' ')
}