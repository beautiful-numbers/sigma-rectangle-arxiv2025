#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bloc A — génération exacte des tables UB_k et des seuils X_def / X_exist
# Sémantique figée :
#   - Famille de premiers : impairs uniquement (2 exclu)
#   - Flux "Bloc A" : q d’abord (q = plus petit premier impair strictement > X),
#                     puis la base 3,5,7,11,13,... en ordre croissant,
#                     en sautant q s’il réapparaît (déduplication).
#   - UB_k = ∏ p/(p−1) sur P_k (premiers distincts)
#   - m_min_S = plus petit m tel que ∏ (1+1/p) ≥ 2 sur le même flux

import json
import math
from fractions import Fraction
from decimal import Decimal
from typing import List, Tuple, Optional, Dict

# ---------- primalité (suffisant pour X ≤ 10^6 et au-delà raisonnable) ----------
def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if (n % 2) == 0:
        return n == 2
    r = int(math.isqrt(n))
    d = 3
    while d <= r:
        if n % d == 0:
            return False
        d += 2
    return True

def next_odd_prime_gt(x: int) -> int:
    """Plus petit premier impair strictement > x."""
    n = x + 1
    if n % 2 == 0:
        n += 1
    if n <= 3:
        n = 3
    while True:
        if is_prime(n):
            return n
        n += 2

def odd_primes_from(start: int = 3):
    """Flux des premiers impairs ≥ start."""
    p = 3 if start <= 3 else (start if start % 2 else start + 1)
    while True:
        if is_prime(p):
            yield p
        p += 2

# ---------- flux Bloc A ----------
def stream_blockA_primes(X: int):
    """
    Suite de premiers DISTINCTS :
      q = plus petit premier impair > X, puis 3,5,7,11,13,... en croissant
      (en sautant q s’il réapparaît).
    """
    q = next_odd_prime_gt(X)
    yield q
    for p in odd_primes_from(3):
        if p == q:
            continue
        yield p

# ---------- produits exacts ----------
def ub_factor(p: int) -> Fraction:
    return Fraction(p, p - 1)   # p/(p-1)

def s_factor_squarefree(p: int) -> Fraction:
    return Fraction(p + 1, p)   # (1 + 1/p)

def product_fraction(ps: List[int], mode: str) -> Fraction:
    acc = Fraction(1, 1)
    if mode == "UB":
        for p in ps:
            acc *= ub_factor(p)
    elif mode == "S":
        for p in ps:
            acc *= s_factor_squarefree(p)
    else:
        raise ValueError("mode inconnu")
    return acc

# ---------- budget en facteurs (digits) ----------
def m_budget(D: int, X: int) -> Tuple[int, List[int]]:
    """
    Max # de facteurs (squarefree) sous digits ≤ D, en empilant le flux Bloc A.
    Critère : ceil(total_log10 + log10 p) ≤ D.
    """
    budget = Decimal(D)
    total = Decimal(0)
    Pk: List[int] = []
    for p in stream_blockA_primes(X):
        nxt = total + Decimal(str(math.log10(p)))
        i = int(nxt)
        ceil_ok = (nxt == Decimal(i)) or (i + 1 <= D)
        if ceil_ok:
            Pk.append(p)
            total = nxt
        else:
            break
    return len(Pk), Pk

def ub_k(X: int, k: int) -> Tuple[Fraction, List[int]]:
    if k <= 0:
        return Fraction(1, 1), []
    ps: List[int] = []
    it = stream_blockA_primes(X)
    while len(ps) < k:
        ps.append(next(it))
    return product_fraction(ps, "UB"), ps

def m_min_S(X: int, cap: int = 200000) -> Tuple[int, List[int], Fraction]:
    acc = Fraction(1, 1)
    plist: List[int] = []
    it = stream_blockA_primes(X)
    for i in range(1, cap + 1):
        p = next(it)
        acc *= s_factor_squarefree(p)
        plist.append(p)
        if acc >= 2:
            return i, plist, acc
    return -1, plist, acc

# ---------- util ----------
def prime_list(a: int, b: int) -> List[int]:
    out: List[int] = []
    for n in range(max(2, a), b + 1):
        if is_prime(n):
            out.append(n)
    return out

def x_def_and_exist(D: int, X_list: List[int]) -> Dict:
    rows = []
    X_def: Optional[int] = None
    X_exist: Optional[int] = None

    for X in X_list:
        k, Pk = m_budget(D, X)
        UB, _ = ub_k(X, k)
        mS, PlS, Sval = m_min_S(X)
        if UB >= 2:
            X_def = X
        if mS != -1 and mS <= k:
            X_exist = X
        rows.append({
            "X": X,
            "k_budget": k,
            "P_k": Pk,
            "UB_k": f"{UB.numerator}/{UB.denominator}",
            "m_min_S": mS,
            "S_at_m_min": f"{Sval.numerator}/{Sval.denominator}",
            "margin_UB": int(2 * UB.denominator - UB.numerator),
        })
    return {"rows": rows, "X_def": X_def, "X_exist": X_exist}

# ---------- main ----------
if __name__ == "__main__":
    # Par défaut : D=6, X premiers de 11 à 200
    D = 6
    Xs = prime_list(11, 200)
    res = x_def_and_exist(D, Xs)
    print(json.dumps(res, ensure_ascii=False, indent=2))
