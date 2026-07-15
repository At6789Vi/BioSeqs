#!/usr/bin/env python3
"""
Python white-box tests for Suffix Array and Suffix Tree.

This implements reference versions of the algorithms for comparison with MoonBit implementation.
Run: python test/python/test_suffix_array_tree.py
"""

def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)

def assert_eq(a, b, msg="Assertion failed"):
    if a != b:
        raise AssertionError(f"{msg}: {a} != {b}")

# ---------------------------------------------------------------------------
# Suffix Array Construction (Prefix-Doubling Algorithm)
# ---------------------------------------------------------------------------

def build_suffix_array(text):
    """Build suffix array using prefix-doubling algorithm"""
    t = text + "$"
    n = len(t)
    sa = list(range(n))
    rank = [ord(c) for c in t]
    k = 1
    while k < n:
        sa.sort(key=lambda i: (rank[i], rank[i + k] if i + k < n else -1))
        new_rank = [0] * n
        r = 0
        new_rank[sa[0]] = r
        for i in range(1, n):
            prev, curr = sa[i - 1], sa[i]
            prev2 = rank[prev + k] if prev + k < n else -1
            curr2 = rank[curr + k] if curr + k < n else -1
            if rank[prev] != rank[curr] or prev2 != curr2:
                r += 1
            new_rank[curr] = r
        rank = new_rank
        k *= 2
    return sa

def get_sorted_suffix(text, sa, idx):
    """Get the suffix at sorted index"""
    pos = sa[idx]
    if pos >= len(text):
        return "$"
    return text[pos:] + "$"

def find_lower(text, sa, pattern):
    """Find lower bound of pattern in suffix array"""
    n = len(sa)
    m = len(pattern)
    left, right = 0, n - 1
    result = -1
    while left <= right:
        mid = (left + right) // 2
        suffix = get_sorted_suffix(text, sa, mid)
        compare = compare_pattern(suffix, pattern, m)
        if compare >= 0:
            if compare == 0:
                result = mid
            right = mid - 1
        else:
            left = mid + 1
    if result == -1 and left < n:
        suffix = get_sorted_suffix(text, sa, left)
        if compare_pattern(suffix, pattern, m) == 0:
            result = left
    return result

def find_upper(text, sa, pattern):
    """Find upper bound of pattern in suffix array"""
    n = len(sa)
    m = len(pattern)
    left, right = 0, n - 1
    result = -1
    while left <= right:
        mid = (left + right) // 2
        suffix = get_sorted_suffix(text, sa, mid)
        compare = compare_pattern(suffix, pattern, m)
        if compare <= 0:
            if compare == 0:
                result = mid
            left = mid + 1
        else:
            right = mid - 1
    return result

def compare_pattern(suffix, pattern, m):
    """Compare suffix with pattern"""
    n = len(suffix)
    for i in range(m):
        if i >= n:
            return -1
        if suffix[i] < pattern[i]:
            return -1
        elif suffix[i] > pattern[i]:
            return 1
    return 0

def sa_count(text, sa, pattern):
    """Count occurrences of pattern"""
    lower = find_lower(text, sa, pattern)
    upper = find_upper(text, sa, pattern)
    if lower == -1 or upper == -1:
        return 0
    return upper - lower + 1

def sa_locate(text, sa, pattern):
    """Locate all occurrences of pattern"""
    lower = find_lower(text, sa, pattern)
    upper = find_upper(text, sa, pattern)
    if lower == -1 or upper == -1:
        return []
    positions = [sa[i] for i in range(lower, upper + 1)]
    return sorted(positions)

# ---------------------------------------------------------------------------
# LCP Array (Longest Common Prefix)
# ---------------------------------------------------------------------------

def build_lcp(text, sa):
    """Build LCP array using Kasai algorithm"""
    t = text + "$"
    n = len(t)
    rank = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    lcp = [0] * n
    k = 0
    for i in range(n - 1):
        j = sa[rank[i] - 1]
        while i + k < n and j + k < n and t[i + k] == t[j + k]:
            k += 1
        lcp[rank[i]] = k
        if k > 0:
            k -= 1
    return lcp

# ---------------------------------------------------------------------------
# Suffix Tree (Naive Implementation)
# ---------------------------------------------------------------------------

def st_contains(text, pattern):
    """Check if pattern is in text"""
    return pattern in text

def st_locate(text, pattern):
    """Find all occurrences of pattern in text"""
    positions = []
    n = len(text)
    m = len(pattern)
    for i in range(n - m + 1):
        if text[i:i + m] == pattern:
            positions.append(i)
    return positions

def st_count(text, pattern):
    """Count occurrences of pattern"""
    return len(st_locate(text, pattern))

def longest_repeated_substring(text):
    """Find longest repeated substring"""
    n = len(text)
    longest = ""
    for i in range(n):
        for j in range(i + 1, n):
            k = 0
            while i + k < n and j + k < n and text[i + k] == text[j + k]:
                k += 1
            if k > len(longest):
                longest = text[i:i + k]
    return longest

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_suffix_array_banana():
    sa = build_suffix_array("banana")
    assert_eq(len(sa), 7)
    assert_eq(get_sorted_suffix("banana", sa, 0), "$")

def test_suffix_array_count():
    sa = build_suffix_array("banana")
    assert_eq(sa_count("banana", sa, "ana"), 2)
    assert_eq(sa_count("banana", sa, "an"), 2)
    assert_eq(sa_count("banana", sa, "b"), 1)
    assert_eq(sa_count("banana", sa, "xyz"), 0)

def test_suffix_array_locate():
    sa = build_suffix_array("banana")
    assert_eq(sa_locate("banana", sa, "ana"), [1, 3])

def test_suffix_array_dna():
    sa = build_suffix_array("ACGTACGT")
    assert_eq(sa_count("ACGTACGT", sa, "ACGT"), 2)
    assert_eq(sa_locate("ACGTACGT", sa, "ACGT"), [0, 4])

def test_suffix_array_empty():
    sa = build_suffix_array("")
    assert_eq(len(sa), 1)

def test_suffix_array_single_char():
    sa = build_suffix_array("a")
    assert_eq(len(sa), 2)
    assert_eq(sa_count("a", sa, "a"), 1)

def test_lcp_banana():
    sa = build_suffix_array("banana")
    lcp = build_lcp("banana", sa)
    assert_eq(len(lcp), 7)
    assert_true(max(lcp) > 0)

def test_lcp_dna():
    sa = build_suffix_array("ACGTACGT")
    lcp = build_lcp("ACGTACGT", sa)
    assert_eq(len(lcp), 9)

def test_suffix_tree_contains():
    assert_true(st_contains("banana", "ana"))
    assert_true(st_contains("banana", "ban"))
    assert_true(not st_contains("banana", "xyz"))

def test_suffix_tree_count():
    assert_eq(st_count("banana", "ana"), 2)
    assert_eq(st_count("banana", "an"), 2)
    assert_eq(st_count("banana", "b"), 1)
    assert_eq(st_count("banana", "xyz"), 0)

def test_suffix_tree_locate():
    assert_eq(st_locate("banana", "ana"), [1, 3])

def test_suffix_tree_longest_repeated():
    assert_eq(longest_repeated_substring("banana"), "ana")
    assert_eq(longest_repeated_substring("aaa"), "aa")

def test_suffix_tree_dna():
    assert_eq(st_count("ACGTACGT", "ACGT"), 2)
    assert_true(st_contains("ACGTACGT", "CGTA"))
    assert_true(not st_contains("ACGTACGT", "CGTT"))

def test_suffix_tree_single_char():
    assert_true(st_contains("a", "a"))
    assert_eq(st_count("a", "a"), 1)

if __name__ == "__main__":
    tests = [
        ("test_suffix_array_banana", test_suffix_array_banana),
        ("test_suffix_array_count", test_suffix_array_count),
        ("test_suffix_array_locate", test_suffix_array_locate),
        ("test_suffix_array_dna", test_suffix_array_dna),
        ("test_suffix_array_empty", test_suffix_array_empty),
        ("test_suffix_array_single_char", test_suffix_array_single_char),
        ("test_lcp_banana", test_lcp_banana),
        ("test_lcp_dna", test_lcp_dna),
        ("test_suffix_tree_contains", test_suffix_tree_contains),
        ("test_suffix_tree_count", test_suffix_tree_count),
        ("test_suffix_tree_locate", test_suffix_tree_locate),
        ("test_suffix_tree_longest_repeated", test_suffix_tree_longest_repeated),
        ("test_suffix_tree_dna", test_suffix_tree_dna),
        ("test_suffix_tree_single_char", test_suffix_tree_single_char),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            test_func()
            print(f"PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL: {name} - {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR: {name} - {type(e).__name__}: {e}")
            failed += 1
    
    print(f"\nTotal: {len(tests)}, Passed: {passed}, Failed: {failed}")
    exit(0 if failed == 0 else 1)
