///|
test "clade_is_terminal" {
  let leaf = Clade::new(name=Some("A"))
  let internal = Clade::new(clades=[Clade::new(name=Some("A"))])
  assert_eq(leaf.is_terminal(), true)
  assert_eq(internal.is_terminal(), false)
}

///|
test "clade_count_terminals" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(clades=[Clade::new(name=Some("C")), Clade::new(name=Some("D"))]),
  ])
  assert_eq(tree.count_terminals(), 4)
}

///|
test "clade_get_terminals" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(clades=[Clade::new(name=Some("C")), Clade::new(name=Some("D"))]),
  ])
  let terminals = tree.get_terminals()
  assert_eq(terminals.length(), 4)
  let names = terminals.map(fn(c) { c.name.unwrap() })
  assert_eq(names.sort(), ["A", "B", "C", "D"].sort())
}

///|
test "clade_get_nonterminals" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(clades=[Clade::new(name=Some("C")), Clade::new(name=Some("D"))]),
  ])
  let nonterminals = tree.get_nonterminals()
  assert_eq(nonterminals.length(), 2)
}

///|
test "clade_find_clades" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(name=Some("Internal"), clades=[
      Clade::new(name=Some("C")),
      Clade::new(name=Some("D")),
    ]),
  ])

  let found = tree.find_clades(name=Some("Internal"))
  assert_eq(found.length(), 1)
  assert_eq(found[0].name, Some("Internal"))

  let terminals = tree.find_clades(terminal=Some(true))
  assert_eq(terminals.length(), 4)
}

///|
test "clade_find_any" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
  ])

  let found = tree.find_any(name=Some("A"))
  assert_eq(found.is_some(), true)
  assert_eq(found.unwrap().name, Some("A"))

  let not_found = tree.find_any(name=Some("Z"))
  assert_eq(not_found.is_none(), true)
}

///|
test "clade_get_path" {
  let tree = Clade::new(name=Some("Root"), clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("Internal"), clades=[
      Clade::new(name=Some("B")),
      Clade::new(name=Some("C")),
    ]),
  ])

  let path = tree.get_path("C")
  assert_eq(path.is_some(), true)
  let names = path.unwrap().map(fn(c) { c.name.unwrap() })
  assert_eq(names.length(), 3)
}

///|
test "clade_depths" {
  let tree = Clade::new(clades=[
    Clade::new(name=Some("A"), branch_length=Some(0.5)),
    Clade::new(name=Some("B"), branch_length=Some(0.3), clades=[
      Clade::new(name=Some("C"), branch_length=Some(0.2)),
    ]),
  ])

  let depths = tree.depths()
  assert_eq(depths.length(), 3)
}

///|
test "clade_is_bifurcating" {
  let bifurcating = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(clades=[Clade::new(name=Some("B")), Clade::new(name=Some("C"))]),
  ])

  let non_bifurcating = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(name=Some("C")),
  ])

  assert_eq(bifurcating.is_bifurcating(), true)
  assert_eq(non_bifurcating.is_bifurcating(), false)
}

///|
test "tree_basic" {
  let root = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
  ])
  let tree = Tree::new(root, rooted=true, name=Some("TestTree"), weight=2.0)

  assert_eq(tree.is_rooted(), true)
  assert_eq(tree.count_terminals(), 2)
  assert_eq(tree.name, Some("TestTree"))
  assert_eq(tree.weight, 2.0)
}

///|
test "tree_get_terminals" {
  let root = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
    Clade::new(clades=[Clade::new(name=Some("C")), Clade::new(name=Some("D"))]),
  ])
  let tree = Tree::new(root)

  let terminals = tree.get_terminals()
  assert_eq(terminals.length(), 4)
}

///|
test "tree_find_clades" {
  let root = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
  ])
  let tree = Tree::new(root)

  let found = tree.find_clades(name=Some("A"))
  assert_eq(found.length(), 1)
}

///|
test "tree_find_any" {
  let root = Clade::new(clades=[
    Clade::new(name=Some("A")),
    Clade::new(name=Some("B")),
  ])
  let tree = Tree::new(root)

  let found = tree.find_any(name=Some("A"))
  assert_eq(found.is_some(), true)
}

///|
test "tree_distance" {
  let root = Clade::new(clades=[
    Clade::new(name=Some("A"), branch_length=Some(0.5)),
    Clade::new(branch_length=Some(0.5), clades=[
      Clade::new(name=Some("B"), branch_length=Some(0.3)),
      Clade::new(name=Some("C"), branch_length=Some(0.3)),
    ]),
  ])
  let tree = Tree::new(root)

  let dist = tree.distance("A", "B")
  assert_eq(dist.is_some(), true)
  assert_eq(dist.unwrap(), 1.3)
}

///|
test "tree_common_ancestor" {
  let root = Clade::new(name=Some("Root"), clades=[
    Clade::new(name=Some("A"), branch_length=Some(0.5)),
    Clade::new(name=Some("Internal"), branch_length=Some(0.5), clades=[
      Clade::new(name=Some("B"), branch_length=Some(0.3)),
      Clade::new(name=Some("C"), branch_length=Some(0.3)),
    ]),
  ])
  let tree = Tree::new(root)

  let mrca = tree.common_ancestor(["B", "C"])
  assert_eq(mrca.is_some(), true)
  assert_eq(mrca.unwrap().name, Some("Internal"))

  let mrca_all = tree.common_ancestor(["A", "B", "C"])
  assert_eq(mrca_all.is_some(), true)
  assert_eq(mrca_all.unwrap().name, Some("Root"))
}

///|
test "newick_parse_simple" {
  let newick = "(A,B);"
  let tree = parse_newick(newick)

  assert_eq(tree.count_terminals(), 2)
  let terminals = tree.get_terminals()
  let names = terminals.map(fn(c) { c.name.unwrap() })
  assert_eq(names.sort(), ["A", "B"].sort())
}

///|
test "newick_parse_nested" {
  let newick = "(((A,B),(C,D)),(E,F,G));"
  let tree = parse_newick(newick)

  assert_eq(tree.count_terminals(), 7)
}

///|
test "newick_parse_with_branch_lengths" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = parse_newick(newick)

  assert_eq(tree.count_terminals(), 4)

  let a = tree.find_any(name=Some("A"))
  assert_eq(a.is_some(), true)
  assert_eq(a.unwrap().branch_length, Some(0.1))

  let c = tree.find_any(name=Some("C"))
  assert_eq(c.is_some(), true)
  assert_eq(c.unwrap().branch_length, Some(0.3))
}

///|
test "newick_parse_no_semicolon" {
  let newick = "(A,B)"
  let tree = parse_newick(newick)

  assert_eq(tree.count_terminals(), 2)
}

///|
test "newick_roundtrip" {
  let original = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = parse_newick(original)
  let result = write_newick(tree)

  let tree2 = parse_newick(result)
  assert_eq(tree2.count_terminals(), 4)
}

///|
test "tree_draw_ascii" {
  let newick = "(A,B,(C,D));"
  let tree = parse_newick(newick)
  let ascii = tree.draw_ascii()

  assert_eq(ascii.contains("A"), true)
  assert_eq(ascii.contains("B"), true)
  assert_eq(ascii.contains("C"), true)
  assert_eq(ascii.contains("D"), true)
}

///|
test "tree_repr" {
  let root = Clade::new(clades=[Clade::new(name=Some("A"))])
  let tree = Tree::new(root, rooted=true, weight=2.5)

  let repr = tree.repr()
  assert_eq(repr.contains("rooted=true"), true)
  assert_eq(repr.contains("weight=2.5"), true)
}
