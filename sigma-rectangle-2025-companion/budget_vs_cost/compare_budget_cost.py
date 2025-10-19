#!/usr/bin/env python3
# compare_budget_cost.py
# Confronte le COÛT D_pivot(y) et le BUDGET D_structure(y, ω) le long d'une expansion simulée.
# - Si --m-requis est fourni, D_pivot = ceil(m_requis * log10(y)).
# - Sinon, calcule une borne sûre D_pivot_lower à partir de la borne constante (rapide).
#
# Expansion = même "esprit" que explore_npi_graph.py (stratégie 'square' par défaut):
#   A: {premier >= y}, B: {q≡1(mod4) >= y, + 2 plus petits premiers >= y, != q}
#   A chaque vague: factorise σ(p^e) et ajoute les nouveaux premiers r >= y (si --respect-y-rude).
# Budget reporté en deux versions à chaque profondeur:
#   - D_struct_min : digits du produit des ω plus petits premiers >= y (borne universelle "worst-case minimale")
#   - D_struct_obs : digits estimés via somme des log10 des premiers réellement rencontrés
#
# Exemples:
#   python compare_budget_cost.py --y 79 --scenario B --max-depth 8 --respect-y-rude --m-requis 912
#   python compare_budget_cost.py --y 101 --scenario A --max-depth 8 --respect-y-rude
#   python compare_budget_cost.py --y 149 --scenario B --max-depth 6 --respect-y-rude --exp-strategy square

import argparse, math, random, sys
from collections import deque

# ---------- utilitaires affichage ----------
def fmt_int(n: int) -> str:
    s = f"{n:,}".replace(",", " ")
    return s

# ---------- primalité / factorisation (MR + Pollard Rho) ----------
_SMALL_PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False
    # Miller–Rabin bases (bonne couverture pratique)
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        ok = False
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                ok = True
                break
        if not ok:
            return False
    return True

def _pollard_rho(n: int) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    if n % 5 == 0:
        return 5
    while True:
        c = random.randrange(1, n-1)
        f = lambda x: (pow(x, 2, n) + c) % n
        x, y, d = random.randrange(2, n-1), random.randrange(2, n-1), 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d

def factorint(n: int, out: dict=None) -> dict:
    if out is None:
        out = {}
    if n == 1:
        return out
    if is_probable_prime(n):
        out[n] = out.get(n, 0) + 1
        return out
    for p in _SMALL_PRIMES:
        if n % p == 0:
            cnt = 0
            while n % p == 0:
                n //= p
                cnt += 1
            out[p] = out.get(p, 0) + cnt
            return factorint(n, out)
    d = _pollard_rho(n)
    factorint(d, out)
    factorint(n // d, out)
    return out

# ---------- premiers / générateurs ----------
def next_prime(n: int) -> int:
    if n <= 2:
        return 2
    p = n + 1 if n % 2 == 0 else n
    while True:
        if is_probable_prime(p):
            return p
        p += 2

def primes_from_y(y: int, m: int):
    # génère les m plus petits premiers >= y
    out = []
    p = y if is_probable_prime(y) else next_prime(y)
    out.append(p)
    while len(out) < m:
        p = next_prime(p + 1)
        out.append(p)
    return out

# ---------- coût: D_pivot ----------
def digits_from_logs(log10_sum: float) -> int:
    return int(log10_sum) + 1 if log10_sum > 0 else 1

def D_pivot_from_m_y(m_requis: int, y: int) -> int:
    return math.ceil(m_requis * math.log10(y))

def lower_bound_m_const(y: int) -> int:
    # borne linéaire m > (y-1)*ln 2
    return int(math.floor((y - 1) * math.log(2))) + 1

def D_pivot_lower_const(y: int) -> int:
    # borne "proof-grade" rapide via m_const
    m0 = lower_bound_m_const(y)
    return math.ceil(m0 * math.log10(y))

# ---------- budget minimal et observé ----------
def D_struct_min(y: int, omega: int) -> int:
    # digits du produit des ω plus petits premiers >= y (via somme des logs)
    if omega <= 0:
        return 1
    plist = primes_from_y(y, omega)
    s = 0.0
    for p in plist:
        s += math.log10(p)
    return digits_from_logs(s)

def D_struct_obs_logsum(primes_set: set) -> int:
    if not primes_set:
        return 1
    s = 0.0
    for p in primes_set:
        s += math.log10(p)
    return digits_from_logs(s)

# ---------- simulateur d'expansion (σ) ----------
def sigma_p_pow_e(p: int, e: int) -> int:
    # σ(p^e) = (p^(e+1)-1)//(p-1)
    return (pow(p, e+1) - 1) // (p - 1)

def init_factors(y: int, scenario: str):
    if scenario.upper() == "A":
        return { next_prime(y) }
    elif scenario.upper() == "B":
        # q ≡ 1 (mod 4) ≥ y, + 2 plus petits >= y (≠ q)
        q = next_prime(y)
        while q % 4 != 1:
            q = next_prime(q + 1)
        p1 = next_prime(y)
        if p1 == q:
            p1 = next_prime(p1 + 1)
        p2 = next_prime(p1 + 1)
        if p2 == q:
            p2 = next_prime(p2 + 1)
        return {q, p1, p2}
    else:
        raise ValueError("scenario must be A or B")

def expand_once(current_frontier: set, y: int, exp_strategy: str, respect_y_rude: bool):
    new_frontier = set()
    newest = None
    for p in list(current_frontier):
        if exp_strategy == "square":
            e = 2 if p % 4 != 1 else 1  # q spécial resté à 1
        else:
            e = 1
        s = sigma_p_pow_e(p, e)
        facs = factorint(s)
        for r in facs.keys():
            if respect_y_rude and r < y:
                continue
            if newest is None or r > newest:
                newest = r
            new_frontier.add(r)
    return new_frontier, newest

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Comparer COÛT (D_pivot) vs BUDGET (D_structure) le long d'une expansion simulée.")
    ap.add_argument("--y", type=int, required=True, help="seuil y-rude (p_min ≥ y)")
    ap.add_argument("--scenario", choices=["A","B"], default="B")
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument("--exp-strategy", choices=["square","p_plus_1"], default="square")
    ap.add_argument("--respect-y-rude", action="store_true")
    ap.add_argument("--m-requis", type=int, help="si fourni: D_pivot = ceil(m_requis*log10(y))")
    ap.add_argument("--print-header", action="store_true")
    args = ap.parse_args()

    y = args.y
    if args.m_requis is not None:
        D_pivot = D_pivot_from_m_y(args.m_requis, y)
        pivot_label = f"D_pivot (m_requis={args.m_requis})"
    else:
        D_pivot = D_pivot_lower_const(y)
        pivot_label = "D_pivot_lower (borne constante)"

    F = init_factors(y, args.scenario)
    frontier = set(F)

    if args.print_header:
        print("\n=== compare_budget_cost ===")
        print(f"y = {y}  | scenario = {args.scenario}  | exp_strategy = {args.exp_strategy}  | max_depth = {args.max_depth}")
        print(f"{pivot_label} = {fmt_int(D_pivot)}\n")

    # table
    print("depth |  omega  |  D_struct_min  |  D_struct_obs  |  newest_prime  |  ratio_min/(ω log10 y)  |  status")
    for d in range(0, args.max_depth + 1):
        omega = len(F)
        Dmin = D_struct_min(y, omega)
        Dobs = D_struct_obs_logsum(F)
        newest = max(F) if F else None
        denom = omega * math.log10(y) if omega > 0 else 1.0
        ratio_min = Dmin / denom if denom > 0 else 0.0
        status = "budget < cost"
        if Dmin >= D_pivot:
            status = "budget ≥ cost"

        np_str = f"{newest}" if newest is not None else "-"
        print(f"{d:5d} | {omega:6d} | {fmt_int(Dmin):14s} | {fmt_int(Dobs):13s} | {np_str:14s} | {ratio_min:20.6f} |  {status}")

        if d == args.max_depth:
            break
        # expansion d -> d+1
        new_frontier, newest_found = expand_once(frontier, y, args.exp_strategy, args.respect_y_rude)
        # union
        to_add = [r for r in new_frontier if r not in F]
        for r in to_add:
            F.add(r)
        frontier = set(to_add)

if __name__ == "__main__":
    # enlever toute limite d'impression d'entiers si Python 3.11+
    if hasattr(sys, "set_int_max_str_digits"):
        try:
            sys.set_int_max_str_digits(0)
        except Exception:
            pass
    main()
