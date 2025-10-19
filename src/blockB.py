#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fichier : src/blockB.py
# Bloc B — Témoins d’abondance (squarefree) en arithmétique exacte.
# Sortie : JSON avec deux sous-blocs :
#   a) familles “avec/sans 3,5,7” (exclusions fixes)
#   b) familles “sans 3,5,7” paramétrées par p_start ∈ {11,13,17,...} (si témoin <= D)

import sys
import json
import math
from fractions import Fraction
from decimal import Decimal
from typing import List, Dict, Set, Tuple, Optional

# --- rendre stdout UTF-8 (évite UnicodeEncodeError sur Windows/redirection) ---
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------- primalité ----------
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

def odd_primes_from(start: int = 3):
    """Flux des premiers impairs >= start."""
    p = 3 if start <= 3 else (start if start % 2 else start + 1)
    while True:
        if is_prime(p):
            yield p
        p += 2

# ---------- produits exacts ----------
def s_factor_squarefree(p: int) -> Fraction:
    return Fraction(p + 1, p)  # (1 + 1/p)

def digits_from_primes(primes: List[int]) -> int:
    """digits(n) = ceil(log10 ∏ p) = ceil(Σ log10 p), p distincts."""
    if not primes:
        return 1
    total = Decimal(0)
    for p in primes:
        total += Decimal(str(math.log10(p)))
    i = int(total)
    return i if total == Decimal(i) else i + 1

def witness_squarefree(D: int, exclusions: Set[int], p_start: int = 3
                      ) -> Tuple[Optional[List[int]], Optional[Fraction]]:
    """
    Construit un témoin squarefree minimal (plus petits premiers admissibles DISTINCTS)
    tel que S ≥ 2 si possible, et vérifie digits <= D.
    Retourne (liste_primes, S_exact) ou (None, None) si aucun témoin <= D.
    """
    S = Fraction(1, 1)
    plist: List[int] = []
    for p in odd_primes_from(p_start):
        if p in exclusions:
            continue
        plist.append(p)
        S *= s_factor_squarefree(p)
        if S >= 2:
            d = digits_from_primes(plist)
            if d <= D:
                return plist, S
            else:
                return None, None

def run_families_a(D: int) -> List[Dict]:
    fams = [
        ("impairs", set()),
        ("sans3", {3}),
        ("sans5", {5}),
        ("sans7", {7}),
        ("sans3_5", {3, 5}),
        ("sans3_7", {3, 7}),
        ("sans5_7", {5, 7}),
    ]
    out = []
    for name, E in fams:
        plist, S = witness_squarefree(D, E, p_start=3)
        if plist is None:
            out.append({
                "family": name,
                "exclusions": sorted(list(E)),
                "witness": None,
                "note": "aucun témoin squarefree <= D",
            })
        else:
            out.append({
                "family": name,
                "exclusions": sorted(list(E)),
                "primes": plist,
                "digits": digits_from_primes(plist),
                "S": f"{S.numerator}/{S.denominator}",
                "margin": f"{(S - 2).numerator}/{(S - 2).denominator}",
            })
    return out

def run_families_b(D: int, pstarts: List[int]) -> List[Dict]:
    out = []
    E = {3, 5, 7}
    for p0 in pstarts:
        plist, S = witness_squarefree(D, E, p_start=p0)
        rec = {
            "p_start": p0,
            "exclusions": sorted(list(E)),
        }
        if plist is None:
            rec["witness"] = None
            rec["note"] = "aucun témoin squarefree <= D"
        else:
            rec.update({
                "primes": plist,
                "digits": digits_from_primes(plist),
                "S": f"{S.numerator}/{S.denominator}",
                "margin": f"{(S - 2).numerator}/{(S - 2).denominator}",
            })
        out.append(rec)
    return out

if __name__ == "__main__":
    # Par défaut : D=6, familles a) et b) (b) avec p_start ∈ {11,13,17,19,23}
    D = 6
    families_a = run_families_a(D)
    pstarts = [11, 13, 17, 19, 23]
    families_b = run_families_b(D, pstarts)
    res = {
        "D": D,
        "families_a": families_a,
        "families_b": families_b,
    }
    print(json.dumps(res, ensure_ascii=False, indent=2))
