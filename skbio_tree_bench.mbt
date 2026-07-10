///|
/// Performance benchmarks for Tree module.

///|
/// Generate a balanced binary tree Newick string with 2^depth tips.
/// Tips are named tip0, tip1, ..., tip(2^depth - 1).
fn gen_balanced_newick(depth : Int) -> String {
  let buf = StringBuilder::new()
  let counter = [0] // mutable counter for tip naming
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
  // internal nodes unnamed (just branch length)
  buf.write_char(':')
  buf.write_string("0.1")
}

///|
/// Generate a caterpillar (ladder) tree Newick string with n tips.
/// E.g., n=4: (((tip0,tip1)int0,tip2)int1,tip3);
fn gen_caterpillar_newick(n : Int) -> String {
  if n < 2 {
    return "tip0:0.1;"
  }
  let buf = StringBuilder::new()
  // n-1 opening parens
  for i = 0; i < n - 1; i = i + 1 {
    buf.write_char('(')
  }
  // First two tips
  buf.write_string("tip0:0.1,tip1:0.1")
  // Close internal nodes and add remaining tips
  for i = 2; i < n; i = i + 1 {
    buf.write_char(')')
    buf.write_string("int")
    buf.write_string((i - 2).to_string())
    buf.write_string(":0.1,tip")
    buf.write_string(i.to_string())
    buf.write_string(":0.1")
  }
  // Final closing paren
  buf.write_char(')')
  buf.write_char(';')
  buf.to_string()
}

///|
test "bench_setup_check" {
  let tree10 = gen_balanced_newick(10)
  let tree = TreeNode::from_newick(tree10)
  assert_eq(tree.count_tips(), 1024)

  let cat = gen_caterpillar_newick(100)
  let cattree = TreeNode::from_newick(cat)
  assert_eq(cattree.count_tips(), 100)
}

///|
test "bench_newick_parse_small" {
  let newick = gen_balanced_newick(8) // 256 tips
  let tree = TreeNode::from_newick(newick)
  assert_eq(tree.count_tips(), 256)
}

///|
test "bench_newick_parse_large" {
  let newick = gen_balanced_newick(12) // 4096 tips
  let tree = TreeNode::from_newick(newick)
  assert_eq(tree.count_tips(), 4096)
}

///|
test "bench_to_newick_small" {
  let newick = gen_balanced_newick(8)
  let tree = TreeNode::from_newick(newick)
  let output = tree.to_newick()
  assert_true(output.has_suffix(";"))
}

///|
test "bench_to_newick_large" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let output = tree.to_newick()
  assert_true(output.has_suffix(";"))
}

///|
test "bench_preorder" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let nodes = tree.preorder()
  assert_eq(nodes.length(), 8191)
}

///|
test "bench_postorder" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let nodes = tree.postorder()
  assert_eq(nodes.length(), 8191)
}

///|
test "bench_levelorder" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let nodes = tree.levelorder()
  assert_eq(nodes.length(), 8191)
}

///|
test "bench_find" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let found = tree.find("tip1000")
  assert_true(found is Some(_))
}

///|
test "bench_find_last" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let found = tree.find("tip4095")
  assert_true(found is Some(_))
}

///|
test "bench_lca" {
  let newick = gen_balanced_newick(10)
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip100")
  names.push("tip200")
  let ancestor = tree.lca(names)
  assert_true(ancestor is Some(_))
}

///|
test "bench_lca_deep" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip0")
  names.push("tip4095")
  let ancestor = tree.lca(names)
  assert_true(ancestor is Some(_))
}

///|
test "bench_distance" {
  let newick = gen_balanced_newick(10)
  let tree = TreeNode::from_newick(newick)
  let dist = tree.distance("tip100", "tip200")
  assert_true(dist is Some(_))
}

///|
test "bench_total_length" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  let total = tree.total_length()
  assert_true(total > 0.0)
}

///|
test "bench_count_tips" {
  let newick = gen_balanced_newick(12)
  let tree = TreeNode::from_newick(newick)
  assert_eq(tree.count_tips(), 4096)
}

///|
test "bench_caterpillar_find" {
  let newick = gen_caterpillar_newick(500)
  let tree = TreeNode::from_newick(newick)
  let found = tree.find("tip498")
  assert_true(found is Some(_))
}

///|
test "bench_caterpillar_lca" {
  let newick = gen_caterpillar_newick(500)
  let tree = TreeNode::from_newick(newick)
  let names = Array::new()
  names.push("tip0")
  names.push("tip498")
  let ancestor = tree.lca(names)
  assert_true(ancestor is Some(_))
}

///|
test "bench_caterpillar_distance" {
  let newick = gen_caterpillar_newick(500)
  let tree = TreeNode::from_newick(newick)
  let dist = tree.distance("tip0", "tip498")
  assert_true(dist is Some(_))
}
