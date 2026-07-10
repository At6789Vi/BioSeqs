#!/usr/bin/env python3
"""
Performance comparison: scikit-bio vs MoonBit Tree module.
Tests the same operations on similar tree structures.
"""

import time
from skbio import TreeNode

def gen_balanced_newick(depth):
    """Generate a balanced binary tree Newick string with 2^depth tips."""
    counter = [0]
    parts = []

    def rec(d):
        if d == 0:
            parts.append(f"tip{counter[0]}:0.1")
            counter[0] += 1
            return
        parts.append("(")
        rec(d - 1)
        parts.append(",")
        rec(d - 1)
        parts.append("):0.1")

    rec(depth)
    parts.append(";")
    return "".join(parts)

def gen_caterpillar_newick(n):
    """Generate a caterpillar (ladder) tree Newick string with n tips."""
    if n < 2:
        return "tip0:0.1;"
    parts = []
    for i in range(n - 1):
        parts.append("(")
    parts.append("tip0:0.1,tip1:0.1")
    for i in range(2, n):
        parts.append(f")int{i-2}:0.1,tip{i}:0.1")
    parts.append(");")
    return "".join(parts)

def time_it(name, func, iterations=1):
    """Time a function and print results."""
    start = time.perf_counter()
    for _ in range(iterations):
        result = func()
    elapsed = time.perf_counter() - start
    print(f"  {name}: {elapsed:.6f}s ({iterations} iterations)")
    return elapsed

def main():
    print("=" * 60)
    print("scikit-bio Tree Performance Benchmark")
    print("=" * 60)
    print()

    # Generate test trees
    newick_1024 = gen_balanced_newick(10)
    newick_4096 = gen_balanced_newick(12)
    cat_500 = gen_caterpillar_newick(500)

    # Parse trees
    tree_1024 = TreeNode.read([newick_1024])
    tree_4096 = TreeNode.read([newick_4096])
    cat_tree = TreeNode.read([cat_500])

    print("--- Balanced Tree (1024 tips) ---")
    time_it("parse", lambda: TreeNode.read([newick_1024]), 10)
    time_it("to_newick", lambda: str(tree_1024), 10)
    time_it("preorder", lambda: list(tree_1024.preorder()), 10)
    time_it("postorder", lambda: list(tree_1024.postorder()), 10)
    time_it("levelorder", lambda: list(tree_1024.levelorder()), 10)
    time_it("find(tip100)", lambda: tree_1024.find("tip100"), 100)
    time_it("find(tip1023)", lambda: tree_1024.find("tip1023"), 100)
    time_it("lca(tip100,tip200)", lambda: tree_1024.lca(["tip100", "tip200"]), 100)
    time_it("distance(tip100,tip200)", lambda: tree_1024.find("tip100").distance(tree_1024.find("tip200")), 100)
    time_it("total_length", lambda: tree_1024.total_length(), 100)
    time_it("count_tips", lambda: len(list(tree_1024.tips())), 100)
    print()

    print("--- Balanced Tree (4096 tips) ---")
    time_it("parse", lambda: TreeNode.read([newick_4096]), 5)
    time_it("to_newick", lambda: str(tree_4096), 5)
    time_it("preorder", lambda: list(tree_4096.preorder()), 5)
    time_it("levelorder", lambda: list(tree_4096.levelorder()), 5)
    time_it("find(tip4095)", lambda: tree_4096.find("tip4095"), 100)
    time_it("lca(tip0,tip4095)", lambda: tree_4096.lca(["tip0", "tip4095"]), 100)
    time_it("distance(tip0,tip4095)", lambda: tree_4096.find("tip0").distance(tree_4096.find("tip4095")), 100)
    time_it("total_length", lambda: tree_4096.total_length(), 50)
    print()

    print("--- Caterpillar Tree (500 tips) ---")
    time_it("parse", lambda: TreeNode.read([cat_500]), 10)
    time_it("find(tip498)", lambda: cat_tree.find("tip498"), 100)
    time_it("lca(tip0,tip498)", lambda: cat_tree.lca(["tip0", "tip498"]), 100)
    time_it("distance(tip0,tip498)", lambda: cat_tree.find("tip0").distance(cat_tree.find("tip498")), 100)
    print()

    print("=" * 60)
    print("Done!")
    print("=" * 60)

if __name__ == "__main__":
    main()
