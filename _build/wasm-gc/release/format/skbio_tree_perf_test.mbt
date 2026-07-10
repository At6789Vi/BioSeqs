///|
/// Performance comparison tests - repeated execution to measure timing.

///|
fn gen_balanced_newick(depth : Int) -> String {
  let buf = StringBuilder::new()
  let counter = [0]
  gen_balanced_rec(buf, depth, counter)
  buf.write_char(';')
  buf.to_string()
}

///|
fn gen_balanced_rec(
  buf : StringBuilder,
  depth : Int,
  counter : Array[Int],
) -> Unit {
  if depth == 0 {
    buf.write_string("tip")
    buf.write_string(counter[0].to_string())
    counter[0] = counter[0] + 1
    buf.write_char(':')
    buf.write_string("0.1")
    return
  }
  buf.write_char('(')
  gen_balanced_rec(buf, depth - 1, counter)
  buf.write_char(',')
  gen_balanced_rec(buf, depth - 1, counter)
  buf.write_char(')')
  buf.write_char(':')
  buf.write_string("0.1")
}

///|
fn gen_caterpillar_newick(n : Int) -> String {
  if n < 2 {
    return "tip0:0.1;"
  }
  let buf = StringBuilder::new()
  for i = 0; i < n - 1; i = i + 1 {
    buf.write_char('(')
  }
  buf.write_string("tip0:0.1,tip1:0.1")
  for i = 2; i < n; i = i + 1 {
    buf.write_char(')')
    buf.write_string("int")
    buf.write_string((i - 2).to_string())
    buf.write_string(":0.1,tip")
    buf.write_string(i.to_string())
    buf.write_string(":0.1")
  }
  buf.write_char(')')
  buf.write_char(';')
  buf.to_string()
}

///|
/// Run LCA 100 times on a 1024-tip balanced tree
test "perf_lca_balanced_1024" {
  let newick = gen_balanced_newick(10) // 1024 tips
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip100")
  names.push("tip200")
  for i = 0; i < 100; i = i + 1 {
    let ancestor = tree.lca(names)
    assert_true(ancestor is Some(_))
  }
}

///|
/// Run LCA 100 times on a 4096-tip balanced tree (deep LCA)
test "perf_lca_balanced_4096_deep" {
  let newick = gen_balanced_newick(12) // 4096 tips
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip0")
  names.push("tip4095")
  for i = 0; i < 100; i = i + 1 {
    let ancestor = tree.lca(names)
    assert_true(ancestor is Some(_))
  }
}

///|
/// Run distance 100 times on a 1024-tip balanced tree
test "perf_distance_balanced_1024" {
  let newick = gen_balanced_newick(10)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 100; i = i + 1 {
    let dist = tree.distance("tip100", "tip200")
    assert_true(dist is Some(_))
  }
}

///|
/// Run distance 100 times on a 4096-tip balanced tree
test "perf_distance_balanced_4096" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 100; i = i + 1 {
    let dist = tree.distance("tip0", "tip4095")
    assert_true(dist is Some(_))
  }
}

///|
/// Run levelorder 50 times on a 4096-tip balanced tree
test "perf_levelorder_4096" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 50; i = i + 1 {
    let nodes = tree.levelorder()
    assert_eq(nodes.length(), 8191)
  }
}

///|
/// Run find 100 times on a 4096-tip balanced tree
test "perf_find_4096" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 100; i = i + 1 {
    let found = tree.find("tip4095")
    assert_true(found is Some(_))
  }
}

///|
/// Run Newick parse 10 times on a 4096-tip balanced tree
test "perf_parse_4096" {
  let newick = gen_balanced_newick(12)
  for i = 0; i < 10; i = i + 1 {
    let tree = TreeNode::from_newick(newick)
    assert_eq(tree.count_tips(), 4096)
  }
}

///|
/// Run to_newick 10 times on a 4096-tip balanced tree
test "perf_to_newick_4096" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 10; i = i + 1 {
    let output = tree.to_newick()
    assert_true(output.has_suffix(";"))
  }
}

///|
/// Run LCA 100 times on a 500-tip caterpillar tree
test "perf_lca_caterpillar_500" {
  let newick = gen_caterpillar_newick(500)
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip0")
  names.push("tip498")
  for i = 0; i < 100; i = i + 1 {
    let ancestor = tree.lca(names)
    assert_true(ancestor is Some(_))
  }
}

///|
/// Run distance 100 times on a 500-tip caterpillar tree
test "perf_distance_caterpillar_500" {
  let newick = gen_caterpillar_newick(500)
  let tree = TreeNode::from_newick(newick)
  for i = 0; i < 100; i = i + 1 {
    let dist = tree.distance("tip0", "tip498")
    assert_true(dist is Some(_))
  }
}
