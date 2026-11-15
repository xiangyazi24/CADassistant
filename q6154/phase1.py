#!/usr/bin/env python3
"""Exact q-series and elimination for Q6154, phase 1.

Discovery script.  It uses SymPy only for rational nullspace/factorisation; all
q-series arithmetic is exact Python integer arithmetic.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import sympy as sp

N = 5000
FIT = 240


def mul_factor(a: list[int], m: int, exponent: int) -> None:
    """In-place multiplication by (1-q^m)^exponent, truncated to len(a)-1."""
    if exponent >= 0:
        for _ in range(exponent):
            for n in range(len(a) - 1, m - 1, -1):
                a[n] -= a[n - m]
    else:
        for _ in range(-exponent):
            for n in range(m, len(a)):
                a[n] += a[n - m]


def eta_quotient_series(residue_exponents: dict[int, int], n: int) -> list[int]:
    """Product over m>=1 with m mod 25 in residue_exponents."""
    a = [0] * (n + 1)
    a[0] = 1
    for m in range(1, n + 1):
        e = residue_exponents.get(m % 25, 0)
        if e:
            mul_factor(a, m, e)
    return a


def shifted(a: list[int], s: int, n: int) -> list[int]:
    out = [0] * (n + 1)
    if s >= 0:
        upto = min(len(a), n + 1 - s)
        out[s:s + upto] = a[:upto]
    else:
        for i in range(-s, min(len(a), n + 1 - s)):
            out[i + s] = a[i]
    return out


def _pack_nonnegative(a: list[int], width: int) -> int:
    z = 0
    for x in reversed(a):
        z = (z << width) + x
    return z


def convolution(a: list[int], b: list[int], n: int) -> list[int]:
    """Exact signed convolution via Kronecker substitution."""
    la = min(len(a), n + 1)
    lb = min(len(b), n + 1)
    while la and a[la - 1] == 0:
        la -= 1
    while lb and b[lb - 1] == 0:
        lb -= 1
    if not la or not lb:
        return [0] * (n + 1)
    aa = a[:la]
    bb = b[:lb]
    # Small cases are faster naively.
    if la * lb <= 20000:
        c = [0] * (n + 1)
        for i, x in enumerate(aa):
            if x:
                lim = min(lb, n + 1 - i)
                for j in range(lim):
                    if bb[j]:
                        c[i + j] += x * bb[j]
        return c
    ap = [max(x, 0) for x in aa]
    an = [max(-x, 0) for x in aa]
    bp = [max(x, 0) for x in bb]
    bn = [max(-x, 0) for x in bb]
    ma = max(max(ap), max(an))
    mb = max(max(bp), max(bn))
    bound = min(la, lb, n + 1) * ma * mb
    width = max(2, bound.bit_length() + 1)
    mask = (1 << width) - 1
    Apos = _pack_nonnegative(ap, width)
    Aneg = _pack_nonnegative(an, width)
    Bpos = _pack_nonnegative(bp, width)
    Bneg = _pack_nonnegative(bn, width)
    vals = [Apos * Bpos, Aneg * Bneg, Apos * Bneg, Aneg * Bpos]
    c = [0] * (n + 1)
    for k in range(n + 1):
        c[k] = (vals[0] & mask) + (vals[1] & mask) - (vals[2] & mask) - (vals[3] & mask)
        vals = [z >> width for z in vals]
    return c


def powers(a: list[int], d: int, n: int) -> list[list[int]]:
    ans = [[1] + [0] * n]
    for _ in range(d):
        ans.append(convolution(ans[-1], a, n))
    return ans


def primitive(v: Iterable[sp.Rational]) -> list[int]:
    vv = list(v)
    den = 1
    for x in vv:
        den = math.lcm(den, int(sp.denom(x)))
    z = [int(x * den) for x in vv]
    g = 0
    for x in z:
        g = math.gcd(g, abs(x))
    z = [x // g for x in z]
    for x in reversed(z):
        if x:
            if x < 0:
                z = [-y for y in z]
            break
    return z


def incremental_rank_mod(rows: list[list[int]], p: int, target: int) -> tuple[int, int]:
    basis: dict[int, list[int]] = {}
    for idx, raw in enumerate(rows):
        v = [x % p for x in raw]
        while True:
            pivot = next((j for j, x in enumerate(v) if x), None)
            if pivot is None:
                break
            if pivot not in basis:
                inv = pow(v[pivot], p - 2, p)
                v = [(x * inv) % p for x in v]
                # Reduce the new pivot out of existing rows to keep a stable RREF.
                for q, w in list(basis.items()):
                    if w[pivot]:
                        t = w[pivot]
                        basis[q] = [(x - t * y) % p for x, y in zip(w, v)]
                basis[pivot] = v
                break
            t = v[pivot]
            w = basis[pivot]
            v = [(x - t * y) % p for x, y in zip(v, w)]
        if len(basis) >= target:
            return idx, len(basis)
    return len(rows) - 1, len(basis)


def digest(a: list[int]) -> str:
    text = ",".join(map(str, a)).encode()
    return hashlib.sha256(text).hexdigest()


def poly_string(coeffs: list[int], monomials: list[tuple[int, int]]) -> str:
    terms: list[str] = []
    for c, (i, j) in zip(coeffs, monomials):
        if not c:
            continue
        atom = ""
        if i:
            atom += "X" if i == 1 else f"X^{i}"
        if j:
            atom += ("*" if atom else "") + ("Y" if j == 1 else f"Y^{j}")
        if not atom:
            atom = "1"
        if not terms:
            if c == -1 and atom != "1":
                terms.append("-" + atom)
            elif c == 1 and atom != "1":
                terms.append(atom)
            else:
                terms.append(f"{c}" if atom == "1" else f"{c}*{atom}")
        else:
            sign = "+" if c > 0 else "-"
            ac = abs(c)
            if ac == 1 and atom != "1":
                terms.append(f" {sign} {atom}")
            else:
                terms.append(f" {sign} {ac}" if atom == "1" else f" {sign} {ac}*{atom}")
    return "".join(terms) if terms else "0"


def main() -> None:
    # E_g contributes at residues g and -g modulo 25.
    def add_e(target: dict[int, int], g: int, e: int) -> None:
        target[g % 25] = target.get(g % 25, 0) + e
        target[(-g) % 25] = target.get((-g) % 25, 0) + e

    eu: dict[int, int] = {}
    for g in (1, 7):
        add_e(eu, g, 1)
    for g in (2, 11):
        add_e(eu, g, -1)

    ev: dict[int, int] = {}
    for g in (1, 7, 3, 4):
        add_e(ev, g, 1)
    for g in (2, 11):
        add_e(ev, g, -2)

    print("building U,V")
    U = shifted(eta_quotient_series(eu, N), 1, N)
    V = shifted(eta_quotient_series(ev, N), 2, N)
    print("building powers")
    Up = powers(U, 7, N)
    Vp = powers(V, 4, N)

    monomials = [(i, j) for j in range(5) for i in range(8)]
    print("building monomial series")
    monos: list[list[int]] = []
    for i, j in monomials:
        if i == 0:
            monos.append(Vp[j])
        elif j == 0:
            monos.append(Up[i])
        else:
            monos.append(convolution(Up[i], Vp[j], N))
        print(f"  monomial X^{i}Y^{j}")

    rows = [[monos[k][n] for k in range(len(monomials))] for n in range(FIT + 1)]
    print("nullspace")
    ns = sp.Matrix(rows).nullspace()
    if len(ns) != 1:
        raise RuntimeError(f"expected one relation, got nullity {len(ns)}")
    coeffs = primitive(ns[0])
    residual = [sum(c * monos[k][n] for k, c in enumerate(coeffs)) for n in range(N + 1)]
    nz = [i for i, x in enumerate(residual) if x]
    if nz:
        raise RuntimeError(f"relation failure first at q^{nz[0]}: {residual[nz[0]]}")

    X, Y = sp.symbols("X Y")
    expr = sum(c * X**i * Y**j for c, (i, j) in zip(coeffs, monomials))
    poly = sp.Poly(expr, X, Y, domain=sp.ZZ)
    factors = sp.factor_list(poly)
    if factors[0] not in (1, -1) or len(factors[1]) != 1 or factors[1][0][1] != 1:
        raise RuntimeError(f"relation reducible: {factors}")

    rows_all = [[monos[k][n] for k in range(len(monomials))] for n in range(FIT + 1)]
    thresholds = []
    for p in (1000000007, 1000000009, 998244353):
        thresholds.append((p,) + incremental_rank_mod(rows_all, p, 39))
    if len({x[1] for x in thresholds}) != 1:
        raise RuntimeError(f"rank thresholds disagree: {thresholds}")
    last_q = thresholds[0][1]
    exact_rank = sp.Matrix(rows_all[:last_q + 1]).rank()
    previous_rank = sp.Matrix(rows_all[:last_q]).rank() if last_q else 0
    if exact_rank != 39 or previous_rank != 38:
        raise RuntimeError((last_q, previous_rank, exact_rank))

    result = {
        "N": N,
        "U_order": next(i for i, x in enumerate(U) if x),
        "V_order": next(i for i, x in enumerate(V) if x),
        "U_coefficients": U,
        "V_coefficients": V,
        "U_sha256_csv": digest(U),
        "V_sha256_csv": digest(V),
        "monomial_order": monomials,
        "relation_coefficients": coeffs,
        "relation": poly_string(coeffs, monomials),
        "relation_sympy": str(poly.as_expr()),
        "degree_X": poly.degree(X),
        "degree_Y": poly.degree(Y),
        "terms": len(poly.terms()),
        "content": int(sp.polys.polytools.terms_gcd(poly.as_expr(), clear=False) == poly.as_expr()),
        "factor_list": str(factors),
        "determining_last_q": last_q,
        "determining_coefficient_count": last_q + 1,
        "rank_before": previous_rank,
        "rank_at": exact_rank,
        "rank_thresholds": thresholds,
        "verified_through_q": N,
        "cleared_common_factor_order": 29,
        "cleared_verified_through_q": N + 29,
        "cusp_orders": {
            "p0": [1, 2],
            "p1": [-1, -3],
            "p2": [2, 1],
            "p3": [1, 4],
            "p4": [-3, -4],
        },
        "pole_degree_U": 4,
        "pole_degree_V": 7,
    }
    out = Path("q6154/phase1.json")
    out.write_text(json.dumps(result, separators=(",", ":")))
    Path("q6154/relation.txt").write_text(result["relation"] + "\n")
    summary = {k: v for k, v in result.items() if k not in ("U_coefficients", "V_coefficients")}
    Path("q6154/phase1-summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("PASS phase1")


if __name__ == "__main__":
    main()
