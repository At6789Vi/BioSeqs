#!/usr/bin/env python3
"""
Python white-box tests for K-means Clustering.

This implements reference versions of the K-means algorithms for comparison with MoonBit implementation.
Run: python test/python/test_kmeans.py
"""

import math

def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)

def assert_eq(a, b, msg="Assertion failed"):
    if a != b:
        raise AssertionError(f"{msg}: {a} != {b}")

# ---------------------------------------------------------------------------
# K-means Clustering Implementation
# ---------------------------------------------------------------------------

def euclidean_distance(a, b):
    """Compute Euclidean distance between two points"""
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

def manhattan_distance(a, b):
    """Compute Manhattan distance between two points"""
    return sum(abs(x - y) for x, y in zip(a, b))

def cosine_distance(a, b):
    """Compute cosine distance between two points"""
    dot = sum(x * y for x, y in zip(a, b))
    norm_a = math.sqrt(sum(x * x for x in a))
    norm_b = math.sqrt(sum(y * y for y in b))
    denom = norm_a * norm_b
    if denom == 0.0:
        return 1.0
    return 1.0 - dot / denom

def assign_clusters(data, centroids):
    """Assign data points to nearest centroid"""
    labels = []
    for point in data:
        min_dist = float('inf')
        min_idx = 0
        for i, centroid in enumerate(centroids):
            dist = euclidean_distance(point, centroid)
            if dist < min_dist:
                min_dist = dist
                min_idx = i
        labels.append(min_idx)
    return labels

def update_centroids(data, labels, k):
    """Update centroids as mean of assigned points"""
    n = len(data)
    if n == 0:
        return []
    
    dim = len(data[0])
    centroids = [[0.0] * dim for _ in range(k)]
    counts = [0] * k
    
    for i, label in enumerate(labels):
        counts[label] += 1
        for j in range(dim):
            centroids[label][j] += data[i][j]
    
    for i in range(k):
        if counts[i] > 0:
            for j in range(dim):
                centroids[i][j] /= counts[i]
    
    return centroids

class KMeans:
    """K-means Clustering implementation"""
    
    def __init__(self, k, max_iter=100, tol=1e-6):
        self.k = k
        self.centroids = []
        self.labels = []
        self.max_iter = max_iter
        self.tol = tol
    
    def fit(self, data):
        """Fit K-means to data"""
        n = len(data)
        if n == 0:
            return self
        
        dim = len(data[0])
        
        centroids = []
        for i in range(self.k):
            idx = i * (n // self.k)
            centroids.append(data[min(idx, n - 1)].copy())
        
        labels = []
        
        for _ in range(self.max_iter):
            new_labels = assign_clusters(data, centroids)
            new_centroids = update_centroids(data, new_labels, self.k)
            
            max_shift = 0.0
            for i in range(self.k):
                shift = euclidean_distance(centroids[i], new_centroids[i])
                if shift > max_shift:
                    max_shift = shift
            
            labels = new_labels
            centroids = new_centroids
            
            if max_shift < self.tol:
                break
        
        self.centroids = centroids
        self.labels = labels
        return self
    
    def predict(self, data):
        """Predict cluster labels for new data"""
        return assign_clusters(data, self.centroids)
    
    def predict_single(self, point):
        """Predict cluster label for a single point"""
        min_dist = float('inf')
        min_idx = 0
        for i, centroid in enumerate(self.centroids):
            dist = euclidean_distance(point, centroid)
            if dist < min_dist:
                min_dist = dist
                min_idx = i
        return min_idx
    
    def get_centroids(self):
        """Get current centroids"""
        return self.centroids
    
    def get_labels(self):
        """Get current labels"""
        return self.labels
    
    def inertia(self, data):
        """Compute sum of squared distances to centroids"""
        total = 0.0
        for i, point in enumerate(data):
            label = self.labels[i]
            dist = euclidean_distance(point, self.centroids[label])
            total += dist * dist
        return total

def silhouette_score(data, labels):
    """Compute silhouette score for clustering"""
    n = len(data)
    if n == 0:
        return 0.0
    
    max_label = max(labels) if labels else 0
    k = max_label + 1
    
    score_sum = 0.0
    
    for i in range(n):
        label = labels[i]
        a = 0.0
        a_count = 0
        b = float('inf')
        
        for j in range(n):
            if i == j:
                continue
            
            dist = euclidean_distance(data[i], data[j])
            
            if labels[j] == label:
                a += dist
                a_count += 1
            else:
                cluster_dist = 0.0
                cluster_count = 0
                for m in range(n):
                    if labels[m] == labels[j]:
                        cluster_dist += euclidean_distance(data[i], data[m])
                        cluster_count += 1
                avg_dist = cluster_dist / cluster_count if cluster_count > 0 else 0.0
                if avg_dist < b:
                    b = avg_dist
        
        a = a / a_count if a_count > 0 else 0.0
        s = 0.0 if max(a, b) == 0.0 else (b - a) / max(a, b)
        score_sum += s
    
    return score_sum / n

def cluster_gene_expression(data, k):
    """Cluster gene expression data"""
    kmeans = KMeans(k, 100, 0.000001)
    return kmeans.fit(data)

def find_optimal_k(data, min_k, max_k):
    """Find optimal k using silhouette score"""
    best_k = min_k
    best_silhouette = -1.0
    
    for k_val in range(min_k, max_k + 1):
        kmeans = KMeans(k_val, 100, 0.000001).fit(data)
        score = silhouette_score(data, kmeans.labels)
        if score > best_silhouette:
            best_silhouette = score
            best_k = k_val
    
    return best_k

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_kmeans_create():
    kmeans = KMeans(3, 100, 0.000001)
    assert_eq(kmeans.k, 3)
    assert_eq(kmeans.max_iter, 100)

def test_euclidean_distance():
    a = [1.0, 2.0]
    b = [4.0, 6.0]
    dist = euclidean_distance(a, b)
    assert_eq(dist, 5.0)

def test_manhattan_distance():
    a = [1.0, 2.0]
    b = [4.0, 6.0]
    dist = manhattan_distance(a, b)
    assert_eq(dist, 7.0)

def test_cosine_distance():
    a = [1.0, 0.0]
    b = [0.0, 1.0]
    dist = cosine_distance(a, b)
    assert_eq(dist, 1.0)

def test_kmeans_fit():
    data = [
        [1.0, 2.0], [2.0, 1.0], [2.0, 3.0],
        [8.0, 7.0], [9.0, 8.0], [7.0, 9.0],
        [4.0, 5.0], [5.0, 4.0], [5.0, 6.0]
    ]
    kmeans = KMeans(3, 100, 0.000001).fit(data)
    assert_eq(len(kmeans.labels), 9)
    assert_eq(len(kmeans.centroids), 3)

def test_kmeans_predict():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0]
    ]
    kmeans = KMeans(2, 100, 0.000001).fit(data)
    new_data = [[1.5, 1.5], [8.5, 8.5]]
    predictions = kmeans.predict(new_data)
    assert_eq(len(predictions), 2)

def test_kmeans_predict_single():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0]
    ]
    kmeans = KMeans(2, 100, 0.000001).fit(data)
    prediction = kmeans.predict_single([1.5, 1.5])
    assert_true(prediction >= 0 and prediction < 2)

def test_kmeans_get_centroids():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0]
    ]
    kmeans = KMeans(2, 100, 0.000001).fit(data)
    centroids = kmeans.get_centroids()
    assert_eq(len(centroids), 2)

def test_kmeans_get_labels():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0]
    ]
    kmeans = KMeans(2, 100, 0.000001).fit(data)
    labels = kmeans.get_labels()
    assert_eq(len(labels), 4)

def test_kmeans_inertia():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0]
    ]
    kmeans = KMeans(2, 100, 0.000001).fit(data)
    inertia = kmeans.inertia(data)
    assert_true(inertia >= 0.0)

def test_kmeans_silhouette_score():
    data = [
        [1.0, 2.0], [2.0, 1.0], [2.0, 3.0],
        [8.0, 7.0], [9.0, 8.0], [7.0, 9.0]
    ]
    labels = [0, 0, 0, 1, 1, 1]
    score = silhouette_score(data, labels)
    assert_true(score > 0.0)

def test_kmeans_cluster_gene_expression():
    data = [
        [1.2, 3.4, 2.1], [2.3, 3.6, 2.3],
        [8.5, 7.2, 6.8], [9.2, 8.1, 7.5],
        [4.2, 5.1, 4.8], [5.3, 6.2, 5.5]
    ]
    kmeans = cluster_gene_expression(data, 3)
    assert_eq(len(kmeans.labels), 6)

def test_kmeans_find_optimal_k():
    data = [
        [1.0, 2.0], [2.0, 1.0],
        [8.0, 7.0], [9.0, 8.0],
        [4.0, 5.0], [5.0, 4.0]
    ]
    k = find_optimal_k(data, 2, 4)
    assert_true(k >= 2 and k <= 4)

def test_kmeans_empty_data():
    data = []
    kmeans = KMeans(3, 100, 0.000001).fit(data)
    assert_eq(len(kmeans.labels), 0)

def test_kmeans_single_point():
    data = [[1.0, 2.0]]
    kmeans = KMeans(1, 100, 0.000001).fit(data)
    assert_eq(len(kmeans.labels), 1)

if __name__ == "__main__":
    tests = [
        ("test_kmeans_create", test_kmeans_create),
        ("test_euclidean_distance", test_euclidean_distance),
        ("test_manhattan_distance", test_manhattan_distance),
        ("test_cosine_distance", test_cosine_distance),
        ("test_kmeans_fit", test_kmeans_fit),
        ("test_kmeans_predict", test_kmeans_predict),
        ("test_kmeans_predict_single", test_kmeans_predict_single),
        ("test_kmeans_get_centroids", test_kmeans_get_centroids),
        ("test_kmeans_get_labels", test_kmeans_get_labels),
        ("test_kmeans_inertia", test_kmeans_inertia),
        ("test_kmeans_silhouette_score", test_kmeans_silhouette_score),
        ("test_kmeans_cluster_gene_expression", test_kmeans_cluster_gene_expression),
        ("test_kmeans_find_optimal_k", test_kmeans_find_optimal_k),
        ("test_kmeans_empty_data", test_kmeans_empty_data),
        ("test_kmeans_single_point", test_kmeans_single_point),
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
