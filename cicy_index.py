"""Numpy/pure-python CICY toolkit: triple intersection numbers, c2.J_r, Euler number, and the index of line bundles.
Configuration: dims = [n_1..n_m] of P^{n_i}; Q = m x K matrix of column degrees q_ij, with sum_j q_ij = n_i + 1 (Calabi-Yau)."""
import itertools, numpy as np
from fractions import Fraction
class CICY:
    def __init__(self, dims, Q):
        self.n = list(dims); self.Q = np.array(Q, int); self.m = len(dims); self.K = self.Q.shape[1]
        assert (self.Q.sum(1) == np.array(dims) + 1).all(), "not Calabi-Yau"
        assert sum(dims) - self.K == 3, "not a threefold"
        self._build()
    # polynomials in J_1..J_m as dict {exponent tuple: Fraction}, truncated by J_i^{n_i+1}=0
    def _mul(self, A, B):
        out = {}
        for ea, ca in A.items():
            for eb, cb in B.items():
                e = tuple(x + y for x, y in zip(ea, eb))
                if any(e[i] > self.n[i] for i in range(self.m)): continue
                out[e] = out.get(e, 0) + ca*cb
        return out
    def _lin(self, coeffs):
        out = {}
        for i, c in enumerate(coeffs):
            if c:
                e = [0]*self.m; e[i] = 1; out[tuple(e)] = Fraction(int(c))
        return out
    def _integrate(self, P):                       # ambient integral: coefficient of prod J_i^{n_i}
        return P.get(tuple(self.n), Fraction(0))
    def _build(self):
        m = self.m; one = {tuple([0]*m): Fraction(1)}
        N = one
        for j in range(self.K): N = self._mul(N, self._lin(self.Q[:, j]))     # normal bundle top class: prod_j (sum_i q_ij J_i)
        self.N = N
        self.d = np.zeros((m, m, m), object)
        for r, s, t in itertools.product(range(m), repeat=3):
            self.d[r, s, t] = self._integrate(self._mul(self._mul(self._mul(self._lin(np.eye(m, dtype=int)[r]), self._lin(np.eye(m, dtype=int)[s])), self._lin(np.eye(m, dtype=int)[t])), N))
        # total Chern class c(TX) = prod_i (1+J_i)^{n_i+1} / prod_j (1 + sum_i q_ij J_i), truncated to degree 3
        def trunc(P, deg): return {e: c for e, c in P.items() if sum(e) <= deg}
        c = one
        for i in range(m):
            f = {tuple([0]*m): Fraction(1)}; f.update(self._lin(np.eye(m, dtype=int)[i]))
            for _ in range(self.n[i] + 1): c = trunc(self._mul(c, f), 3)
        for j in range(self.K):
            x = self._lin(self.Q[:, j]); inv = dict(one); term = dict(one)
            for p in range(1, 4):
                term = trunc(self._mul(term, x), 3); sign = (-1)**p
                for e, v in term.items(): inv[e] = inv.get(e, 0) + sign*v
            c = trunc(self._mul(c, inv), 3)
        self.c2J = np.array([self._integrate(self._mul(self._mul({e: v for e, v in c.items() if sum(e) == 2}, self._lin(np.eye(m, dtype=int)[r])), N)) for r in range(m)], object)
        self.euler = self._integrate(self._mul({e: v for e, v in c.items() if sum(e) == 3}, N))
    def index(self, k):
        k = np.array(k, int)
        d3 = sum(self.d[r, s, t]*k[r]*k[s]*k[t] for r, s, t in itertools.product(range(self.m), repeat=3))
        return d3/6 + sum(self.c2J[r]*k[r] for r in range(self.m))/12
if __name__ == "__main__":
    tests = {"quintic [P4|5]": ([4], [[5]], (-200, 50)),
             "bicubic [P2 P2|3 3]": ([2, 2], [[3], [3]], (-162, None)),
             "[P1|2][P3|4]": ([1, 3], [[2], [4]], (-168, None)),
             "tetraquadric (P1)^4|(2,2,2,2)": ([1, 1, 1, 1], [[2], [2], [2], [2]], (-128, None)),
             }
    for name, t in tests.items():
        if t is None: continue
        dims, Q, ex = t; X = CICY(dims, Q)
        print(f"{name:34s} euler={X.euler}  (expected {ex[0]})  c2.J={list(X.c2J)}  d111={X.d[0,0,0]}")
    q = CICY([4], [[5]])
    print("quintic index O(1) =", q.index([1]), "(expect h0=5);  O(-1) =", q.index([-1]), " O(2)=", q.index([2]), "(expect 15 = h0)")
