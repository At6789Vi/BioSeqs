///|
/// Tests for Bio.PDB module.

///|
test "vector3_distance" {
  let v1 = Vector3::new(0.0, 0.0, 0.0)
  let v2 = Vector3::new(3.0, 4.0, 0.0)
  assert_eq(v1.distance(v2), 5.0)
}

///|
test "atom_creation" {
  let coord = Vector3::new(1.0, 2.0, 3.0)
  let atom = Atom::new(name="CA", coord=coord, resname="ALA", chainid='A', resseq=1)
  
  assert_eq(atom.name, "CA")
  assert_eq(atom.resname, "ALA")
  assert_eq(atom.chainid, 'A')
  assert_eq(atom.resseq, 1)
  assert_eq(atom.coord.x, 1.0)
}

///|
test "atom_distance" {
  let atom1 = Atom::new(name="CA", coord=Vector3::new(0.0, 0.0, 0.0))
  let atom2 = Atom::new(name="CB", coord=Vector3::new(1.5, 0.0, 0.0))
  
  assert_eq(atom1.distance(atom2), 1.5)
}

///|
test "residue_creation" {
  let atoms = [
    Atom::new(name="N", coord=Vector3::new(0.0, 0.0, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(1.5, 0.0, 0.0)),
    Atom::new(name="C", coord=Vector3::new(3.0, 0.0, 0.0)),
  ]
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms)
  
  assert_eq(residue.resname, "ALA")
  assert_eq(residue.get_num_atoms(), 3)
}

///|
test "residue_find_atom" {
  let atoms = [
    Atom::new(name="N", coord=Vector3::new(0.0, 0.0, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(1.5, 0.0, 0.0)),
  ]
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms)
  
  let ca = residue.find_atom("CA")
  assert_eq(ca is Some(_), true)
  assert_eq(ca.unwrap().name, "CA")
  
  let cb = residue.find_atom("CB")
  assert_eq(cb is None, true)
}

///|
test "residue_get_ca" {
  let atoms = [
    Atom::new(name="N", coord=Vector3::new(0.0, 0.0, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(1.5, 0.0, 0.0)),
  ]
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms)
  
  let ca = residue.get_ca()
  assert_eq(ca is Some(_), true)
}

///|
test "residue_center_of_mass" {
  let atoms = [
    Atom::new(name="N", coord=Vector3::new(0.0, 0.0, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(2.0, 0.0, 0.0)),
    Atom::new(name="C", coord=Vector3::new(4.0, 0.0, 0.0)),
  ]
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms)
  
  let com = residue.center_of_mass()
  assert_eq(com is Some(_), true)
  assert_eq(com.unwrap().x, 2.0)
}

///|
test "chain_creation" {
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1)
  let chain = Chain::new(id='A', residues=[residue])
  
  assert_eq(chain.get_id(), 'A')
  assert_eq(chain.get_num_residues(), 1)
}

///|
test "chain_find_residue" {
  let r1 = Residue::new(resname="ALA", chainid='A', resseq=1)
  let r2 = Residue::new(resname="SER", chainid='A', resseq=2)
  let chain = Chain::new(id='A', residues=[r1, r2])
  
  let found = chain.find_residue(2)
  assert_eq(found is Some(_), true)
  assert_eq(found.unwrap().resname, "SER")
}

///|
test "model_creation" {
  let chain = Chain::new(id='A')
  let model = Model::new(id=1, chains=[chain])
  
  assert_eq(model.get_id(), 1)
  assert_eq(model.get_num_chains(), 1)
}

///|
test "model_find_chain" {
  let chain_a = Chain::new(id='A')
  let chain_b = Chain::new(id='B')
  let model = Model::new(id=1, chains=[chain_a, chain_b])
  
  let found = model.find_chain('B')
  assert_eq(found is Some(_), true)
}

///|
test "structure_creation" {
  let model = Model::new(id=1)
  let structure = Structure::new(id="test", models=[model])
  
  assert_eq(structure.get_id(), "test")
  assert_eq(structure.get_num_models(), 1)
}

///|
test "pdb_parse_simple" {
  let pdb_content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\nATOM      3  C   ALA A   1      3.000   0.000   0.000  1.00  0.00           C\nEND\n"
  
  let structure = parse_pdb(pdb_content)
  
  assert_eq(structure.get_num_models(), 1)
  assert_eq(structure.get_num_chains(), 1)
  assert_eq(structure.get_num_residues(), 1)
  assert_eq(structure.get_num_atoms(), 3)
}

///|
test "pdb_parse_multiple_residues" {
  let pdb_content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\nATOM      3  N   SER A   2      3.000   0.000   0.000  1.00  0.00           N\nATOM      4  CA  SER A   2      4.500   0.000   0.000  1.00  0.00           C\nEND\n"
  
  let structure = parse_pdb(pdb_content)
  
  assert_eq(structure.get_num_residues(), 2)
  assert_eq(structure.get_num_atoms(), 4)
}

///|
test "pdb_parse_multiple_chains" {
  let pdb_content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\nATOM      3  N   SER B   1      5.000   0.000   0.000  1.00  0.00           N\nATOM      4  CA  SER B   1      6.500   0.000   0.000  1.00  0.00           C\nEND\n"
  
  let structure = parse_pdb(pdb_content)
  
  assert_eq(structure.get_num_chains(), 2)
}

///|
test "pdb_parse_with_model" {
  let pdb_content = "MODEL     1\nATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\nENDMDL\nEND\n"
  
  let structure = parse_pdb(pdb_content)
  
  assert_eq(structure.get_num_models(), 1)
}

///|
test "pdb_write" {
  let atom = Atom::new(name="CA", coord=Vector3::new(1.0, 2.0, 3.0), resname="ALA", chainid='A', resseq=1)
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[atom])
  let chain = Chain::new(id='A', residues=[residue])
  let model = Model::new(id=1, chains=[chain])
  let structure = Structure::new(id="test", models=[model])
  
  let pdb_str = write_pdb(structure)
  
  assert_eq(pdb_str.contains("ATOM"), true)
  assert_eq(pdb_str.contains("ALA"), true)
}

///|
test "pdb_roundtrip" {
  let pdb_content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\nEND\n"
  
  let structure = parse_pdb(pdb_content)
  let written = write_pdb(structure)
  let re_parsed = parse_pdb(written)
  
  assert_eq(re_parsed.get_num_atoms(), 2)
}

///|
test "structure_distance" {
  let atom1 = Atom::new(name="CA", coord=Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1)
  let atom2 = Atom::new(name="CB", coord=Vector3::new(3.0, 4.0, 0.0), resname="ALA", chainid='A', resseq=1)
  let residue = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[atom1, atom2])
  let chain = Chain::new(id='A', residues=[residue])
  let model = Model::new(id=1, chains=[chain])
  let structure = Structure::new(id="test", models=[model])
  
  let dist = structure.distance("CA", "CB")
  assert_eq(dist is Some(_), true)
  assert_eq(dist.unwrap(), 5.0)
}