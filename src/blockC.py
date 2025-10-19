#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bloc C — Branch-and-Bound minimal (coupes/validations exactes) pour D modeste.
# Objectif : produire des exemples représentatifs (pas une exploration exhaustive lourde).
# - Arithmétique 100% rationnelle pour S, UB_headroom, UB_new, Smin_rest.
# - Budget par digits(n) = ceil(sum e_i log10 p_i) ≤ D.
# - Primes impairs uniquement.

import json
import math
from fractions import Fraction
from decimal import Decimal
from typing import List, Dict, Tuple

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

def prime_list(limit: int) -> List[int]:
    out = []
    for n in range(3, limit + 1, 2):
        if is_prime(n):
            out.append(n)
    return out

# ---------- arithmétique locale ----------
def sigma_over_p_power(p: int, e: int) -> Fraction:
    # sigma(p^e)/p^e = (p^(e+1)-1)/(p^e*(p-1))
    num = pow(p, e + 1) - 1
    den = pow(p, e) * (p - 1)
    return Fraction(num, den)

def ub_factor(p: int) -> Fraction:
    return Fraction(p, p - 1)  # p/(p-1)

def s_factor_squarefree(p: int) -> Fraction:
    return Fraction(p + 1, p)  # (1 + 1/p)

# ---------- util digits ----------
def ceil_log10_mul(current_log10: Decimal, p: int, e_add: int = 1) -> int:
    nxt = current_log10 + Decimal(str(e_add * math.log10(p)))
    i = int(nxt)
    return i if nxt == Decimal(i) else i + 1

# ---------- BnB minimal ----------
class BnB:
    def __init__(self, D: int, primes: List[int], max_nodes: int = 50000):
        self.D = D
        self.primes = primes
        self.max_nodes = max_nodes
        self.nodes = 0
        self.cuts: List[Dict] = []
        self.validations: List[Dict] = []
        self.equalities: List[Dict] = []

    def explore(self):
        self._dfs(
            idx=0,
            log10_n=Decimal(0),
            S_part=Fraction(1, 1),
            factors=[],
        )

    def _dfs(self, idx: int, log10_n: Decimal, S_part: Fraction, factors: List[Tuple[int,int]]):
        if self.nodes >= self.max_nodes:
            return
        self.nodes += 1

        # Construire UB_headroom (sur premiers déjà présents) et UB_new (sur nouveaux)
        UB_headroom = Fraction(1, 1)
        for (p, e) in factors:
            UB_headroom *= ub_factor(p) / sigma_over_p_power(p, e)

        # UB_new: p/(p-1) sur nouveaux premiers tant que digits ≤ D
        UB_new = Fraction(1, 1)
        log_tmp = Decimal(log10_n)  # copie
        used = {p for (p, _) in factors}
        for j in range(idx, len(self.primes)):
            p = self.primes[j]
            if p in used:
                continue
            # test digits si on ajoute p^1
            if ceil_log10_mul(log_tmp, p, 1) > self.D:
                break
            UB_new *= ub_factor(p)
            log_tmp += Decimal(str(math.log10(p)))

        # Coupe sûre ?
        if S_part * UB_headroom * UB_new < 2:
            self.cuts.append({
                "factors": factors.copy(),
                "S_part": f"{S_part.numerator}/{S_part.denominator}",
                "UB_headroom": f"{UB_headroom.numerator}/{UB_headroom.denominator}",
                "UB_new": f"{UB_new.numerator}/{UB_new.denominator}",
                "decision": "cut",
            })
            return

        # Validation (Smin_rest) : (1+1/p) sur nouveaux, tant que digits ≤ D
        Smin_rest = Fraction(1, 1)
        log_tmp2 = Decimal(log10_n)
        used2 = {p for (p, _) in factors}
        planned: List[int] = []
        for j in range(idx, len(self.primes)):
            p = self.primes[j]
            if p in used2:
                continue
            if ceil_log10_mul(log_tmp2, p, 1) > self.D:
                break
            Smin_rest *= s_factor_squarefree(p)
            log_tmp2 += Decimal(str(math.log10(p)))
            planned.append(p)

        if S_part * Smin_rest > 2:
            # Construire un témoin : factors + planned (exposants 1)
            witness = factors.copy() + [(p, 1) for p in planned]
            S_full = Fraction(1, 1)
            for (p, e) in witness:
                S_full *= sigma_over_p_power(p, e)
            self.validations.append({
                "witness": witness,
                "S": f"{S_full.numerator}/{S_full.denominator}",
                "status": "abundant",
            })
            return

        # Feuille ? (plus de primes ou budget bloqué)
        if idx >= len(self.primes):
            S_full = S_part
            status = "perfect" if S_full == 2 else ("abundant" if S_full > 2 else "deficient")
            if status == "perfect":
                self.equalities.append({
                    "witness": factors.copy(),
                    "S": f"{S_full.numerator}/{S_full.denominator}",
                })
            return

        # Branche 1 : exclure p_idx
        self._dfs(idx + 1, log10_n, S_part, factors)

        # Branche 2 : inclure p_idx avec e = 1,2,... tant que digits ≤ D
        p = self.primes[idx]
        e = 1
        log_local = Decimal(log10_n)
        S_local = Fraction(S_part)
        while True:
            # test digits
            if ceil_log10_mul(log_local, p, 1) > self.D:
                break
            log_local += Decimal(str(math.log10(p)))
            S_local *= sigma_over_p_power(p, e)
            self._dfs(idx + 1, log_local, S_local, factors + [(p, e)])
            e += 1
            if self.nodes >= self.max_nodes:
                break

if __name__ == "__main__":
    # Par défaut : D=6, premiers de travail jusqu’à 97 (impairs)
    D = 6
    primes = prime_list(97)  # [3,5,7,11,13,17,19,23,29,31, ..., 97]

    bnb = BnB(D=D, primes=primes, max_nodes=30000)
    bnb.explore()

    res = {
        "D": D,
        "nodes_visited": bnb.nodes,
        "cuts": bnb.cuts[:20],          # échantillon pour ne pas inonder
        "validations": bnb.validations, # témoins trouvés
        "equalities": bnb.equalities,   # cas S=2 si jamais
    }
    print(json.dumps(res, ensure_ascii=False, indent=2))
