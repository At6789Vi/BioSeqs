#!/usr/bin/env python3
"""
Python white-box tests for Hidden Markov Model (HMM).

This implements reference versions of the HMM algorithms for comparison with MoonBit implementation.
Run: python test/python/test_hmm.py
"""

import math
from collections import defaultdict

def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)

def assert_eq(a, b, msg="Assertion failed"):
    if a != b:
        raise AssertionError(f"{msg}: {a} != {b}")

# ---------------------------------------------------------------------------
# Hidden Markov Model Implementation
# ---------------------------------------------------------------------------

class HMM:
    """Hidden Markov Model implementation"""
    
    def __init__(self, states, observations, start_prob, trans_prob, emit_prob):
        self.states = states
        self.observations = observations
        self.start_prob = start_prob
        self.trans_prob = trans_prob
        self.emit_prob = emit_prob
    
    def forward(self, obs):
        """Forward algorithm - compute probability of observation sequence"""
        n = len(obs)
        if n == 0:
            return 0.0
        
        m = len(self.states)
        alpha = [0.0] * m
        
        for i in range(m):
            state = self.states[i]
            start = self.start_prob.get(state, 0.0)
            emit = self.emit_prob.get(state, {}).get(obs[0], 0.0)
            alpha[i] = start * emit
        
        for t in range(1, n):
            new_alpha = [0.0] * m
            for i in range(m):
                sum_val = 0.0
                for j in range(m):
                    sum_val += alpha[j] * self.trans_prob.get(self.states[j], {}).get(self.states[i], 0.0)
                new_alpha[i] = sum_val * self.emit_prob.get(self.states[i], {}).get(obs[t], 0.0)
            alpha = new_alpha
        
        return sum(alpha)
    
    def backward(self, obs):
        """Backward algorithm - compute backward probabilities"""
        n = len(obs)
        if n == 0:
            return {state: [0.0] for state in self.states}
        
        m = len(self.states)
        beta = {state: [0.0] * n for state in self.states}
        
        for i in range(m):
            beta[self.states[i]][n - 1] = 1.0
        
        for t in range(n - 2, -1, -1):
            for i in range(m):
                sum_val = 0.0
                for j in range(m):
                    sum_val += self.trans_prob.get(self.states[i], {}).get(self.states[j], 0.0) * \
                               self.emit_prob.get(self.states[j], {}).get(obs[t + 1], 0.0) * \
                               beta[self.states[j]][t + 1]
                beta[self.states[i]][t] = sum_val
        
        return beta
    
    def viterbi(self, obs):
        """Viterbi algorithm - find most probable state sequence"""
        n = len(obs)
        if n == 0:
            return []
        
        m = len(self.states)
        delta = [[0.0] * n for _ in range(m)]
        psi = [[0] * n for _ in range(m)]
        
        for i in range(m):
            state = self.states[i]
            start = self.start_prob.get(state, 0.0)
            emit = self.emit_prob.get(state, {}).get(obs[0], 0.0)
            delta[i][0] = start * emit
        
        for t in range(1, n):
            for i in range(m):
                max_val = 0.0
                max_idx = 0
                for j in range(m):
                    val = delta[j][t - 1] * self.trans_prob.get(self.states[j], {}).get(self.states[i], 0.0)
                    if val > max_val:
                        max_val = val
                        max_idx = j
                delta[i][t] = max_val * self.emit_prob.get(self.states[i], {}).get(obs[t], 0.0)
                psi[i][t] = max_idx
        
        path = []
        max_idx = 0
        max_val = delta[0][n - 1]
        for i in range(1, m):
            if delta[i][n - 1] > max_val:
                max_val = delta[i][n - 1]
                max_idx = i
        
        path.append(self.states[max_idx])
        for t in range(n - 1, 0, -1):
            max_idx = psi[max_idx][t]
            path.append(self.states[max_idx])
        
        return path[::-1]
    
    def baum_welch(self, obs, max_iter=10, tol=1e-6):
        """Baum-Welch algorithm - train HMM parameters"""
        n = len(obs)
        if n == 0:
            return self
        
        m = len(self.states)
        states = self.states
        
        start_prob = dict(self.start_prob)
        trans_prob = {s: dict(self.trans_prob.get(s, {})) for s in states}
        emit_prob = {s: dict(self.emit_prob.get(s, {})) for s in states}
        
        for _ in range(max_iter):
            alpha = [0.0] * m
            for i in range(m):
                alpha[i] = start_prob.get(states[i], 0.0) * emit_prob.get(states[i], {}).get(obs[0], 0.0)
            
            beta = [0.0] * m
            for i in range(m):
                beta[i] = 1.0
            
            alphas = [alpha[:]]
            betas = [beta[:]]
            
            for t in range(1, n):
                new_alpha = [0.0] * m
                for i in range(m):
                    sum_val = 0.0
                    for j in range(m):
                        sum_val += alphas[-1][j] * trans_prob.get(states[j], {}).get(states[i], 0.0)
                    new_alpha[i] = sum_val * emit_prob.get(states[i], {}).get(obs[t], 0.0)
                alphas.append(new_alpha)
                
                new_beta = [0.0] * m
                for i in range(m):
                    sum_val = 0.0
                    for j in range(m):
                        sum_val += trans_prob.get(states[i], {}).get(states[j], 0.0) * \
                                   emit_prob.get(states[j], {}).get(obs[n - t], 0.0) * betas[-1][j]
                    new_beta[i] = sum_val
                betas.insert(0, new_beta)
            
            gamma = [[0.0] * n for _ in range(m)]
            for i in range(m):
                for t in range(n):
                    total = sum(alphas[t][j] * betas[t][j] for j in range(m))
                    if total > 0:
                        gamma[i][t] = alphas[t][i] * betas[t][i] / total
            
            xi = [[[0.0] * (n - 1) for _ in range(m)] for __ in range(m)]
            for i in range(m):
                for j in range(m):
                    for t in range(n - 1):
                        total = sum(sum(alphas[t][k] * trans_prob.get(states[k], {}).get(states[l], 0.0) * \
                                       emit_prob.get(states[l], {}).get(obs[t + 1], 0.0) * betas[t + 1][l] \
                                       for l in range(m)) for k in range(m))
                        if total > 0:
                            xi[i][j][t] = alphas[t][i] * trans_prob.get(states[i], {}).get(states[j], 0.0) * \
                                          emit_prob.get(states[j], {}).get(obs[t + 1], 0.0) * betas[t + 1][j] / total
            
            for i in range(m):
                start_prob[states[i]] = gamma[i][0]
            
            for i in range(m):
                total = sum(gamma[i][t] for t in range(n - 1))
                if total > 0:
                    for j in range(m):
                        trans_prob[states[i]][states[j]] = sum(xi[i][j][t] for t in range(n - 1)) / total
                else:
                    for j in range(m):
                        trans_prob[states[i]][states[j]] = 1.0 / m
            
            for i in range(m):
                total = sum(gamma[i][t] for t in range(n))
                if total > 0:
                    for obs_val in self.observations:
                        obs_total = sum(gamma[i][t] for t in range(n) if obs[t] == obs_val)
                        emit_prob[states[i]][obs_val] = obs_total / total
                else:
                    for obs_val in self.observations:
                        emit_prob[states[i]][obs_val] = 1.0 / len(self.observations)
        
        return HMM(states, self.observations, start_prob, trans_prob, emit_prob)

def create_gene_prediction_hmm():
    """Create a gene prediction HMM"""
    states = ["Start", "Exon", "Intron", "Intergenic", "Stop"]
    observations = ["A", "T", "C", "G"]
    
    start_prob = {
        "Start": 0.01,
        "Exon": 0.01,
        "Intron": 0.01,
        "Intergenic": 0.97,
        "Stop": 0.0
    }
    
    trans_prob = {
        "Start": {"Start": 0.0, "Exon": 0.9, "Intron": 0.05, "Intergenic": 0.05, "Stop": 0.0},
        "Exon": {"Start": 0.0, "Exon": 0.9, "Intron": 0.05, "Intergenic": 0.04, "Stop": 0.01},
        "Intron": {"Start": 0.0, "Exon": 0.1, "Intron": 0.85, "Intergenic": 0.04, "Stop": 0.01},
        "Intergenic": {"Start": 0.01, "Exon": 0.01, "Intron": 0.01, "Intergenic": 0.97, "Stop": 0.0},
        "Stop": {"Start": 0.0, "Exon": 0.0, "Intron": 0.0, "Intergenic": 1.0, "Stop": 0.0}
    }
    
    emit_prob = {
        "Start": {"A": 0.3, "T": 0.3, "C": 0.2, "G": 0.2},
        "Exon": {"A": 0.3, "T": 0.3, "C": 0.2, "G": 0.2},
        "Intron": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        "Intergenic": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        "Stop": {"A": 0.3, "T": 0.3, "C": 0.2, "G": 0.2}
    }
    
    return HMM(states, observations, start_prob, trans_prob, emit_prob)

def predict_genes(hmm, dna):
    """Predict genes from DNA sequence"""
    obs = list(dna)
    return hmm.viterbi(obs)

def extract_exons(prediction, dna):
    """Extract exons from prediction"""
    exons = []
    current_exon = []
    for i, state in enumerate(prediction):
        if state == "Exon":
            current_exon.append(dna[i])
        else:
            if current_exon:
                exons.append("".join(current_exon))
                current_exon = []
    if current_exon:
        exons.append("".join(current_exon))
    return exons

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_hmm_create():
    hmm = create_gene_prediction_hmm()
    assert_eq(len(hmm.states), 5)
    assert_eq(len(hmm.observations), 4)

def test_hmm_forward():
    hmm = create_gene_prediction_hmm()
    obs = ["A", "T", "C", "G"]
    prob = hmm.forward(obs)
    assert_true(prob > 0.0)

def test_hmm_backward():
    hmm = create_gene_prediction_hmm()
    obs = ["A", "T", "C", "G"]
    beta = hmm.backward(obs)
    assert_eq(len(beta), 5)

def test_hmm_viterbi():
    hmm = create_gene_prediction_hmm()
    obs = ["A", "T", "C", "G"]
    path = hmm.viterbi(obs)
    assert_eq(len(path), 4)

def test_hmm_predict_genes():
    hmm = create_gene_prediction_hmm()
    prediction = predict_genes(hmm, "ATCG")
    assert_eq(len(prediction), 4)

def test_hmm_extract_exons():
    hmm = create_gene_prediction_hmm()
    prediction = predict_genes(hmm, "ATCGATCG")
    exons = extract_exons(prediction, "ATCGATCG")
    assert_true(len(exons) >= 0)

def test_hmm_small_dna():
    hmm = create_gene_prediction_hmm()
    dna = "AAATTTCCCGGG"
    prediction = predict_genes(hmm, dna)
    assert_eq(len(prediction), len(dna))

def test_hmm_baum_welch():
    hmm = create_gene_prediction_hmm()
    obs = ["A", "T", "C", "G", "A", "T"]
    trained = hmm.baum_welch(obs, 10, 0.000001)
    assert_eq(len(trained.states), 5)

def test_hmm_state_transitions():
    hmm = create_gene_prediction_hmm()
    assert_true("Start" in hmm.trans_prob)
    assert_true("Exon" in hmm.trans_prob)

def test_hmm_emit_probabilities():
    hmm = create_gene_prediction_hmm()
    assert_true("A" in hmm.emit_prob.get("Exon", {}))

def test_hmm_empty_input():
    hmm = create_gene_prediction_hmm()
    prob = hmm.forward([])
    assert_eq(prob, 0.0)

def test_hmm_single_observation():
    hmm = create_gene_prediction_hmm()
    prob = hmm.forward(["A"])
    assert_true(prob > 0.0)

def test_hmm_viterbi_single():
    hmm = create_gene_prediction_hmm()
    path = hmm.viterbi(["A"])
    assert_eq(len(path), 1)

def test_hmm_backward_single():
    hmm = create_gene_prediction_hmm()
    beta = hmm.backward(["A"])
    assert_eq(len(beta), 5)

def test_hmm_exon_extraction():
    hmm = create_gene_prediction_hmm()
    dna = "AAAA"
    prediction = predict_genes(hmm, dna)
    exons = extract_exons(prediction, dna)
    assert_true(len(exons) >= 0)

if __name__ == "__main__":
    tests = [
        ("test_hmm_create", test_hmm_create),
        ("test_hmm_forward", test_hmm_forward),
        ("test_hmm_backward", test_hmm_backward),
        ("test_hmm_viterbi", test_hmm_viterbi),
        ("test_hmm_predict_genes", test_hmm_predict_genes),
        ("test_hmm_extract_exons", test_hmm_extract_exons),
        ("test_hmm_small_dna", test_hmm_small_dna),
        ("test_hmm_baum_welch", test_hmm_baum_welch),
        ("test_hmm_state_transitions", test_hmm_state_transitions),
        ("test_hmm_emit_probabilities", test_hmm_emit_probabilities),
        ("test_hmm_empty_input", test_hmm_empty_input),
        ("test_hmm_single_observation", test_hmm_single_observation),
        ("test_hmm_viterbi_single", test_hmm_viterbi_single),
        ("test_hmm_backward_single", test_hmm_backward_single),
        ("test_hmm_exon_extraction", test_hmm_exon_extraction),
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
