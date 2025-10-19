#!/usr/bin/env python3
# explore_npi_graph.py
# Étape 5a : exploration empirique de l'expansion forcée via r | σ(p^e)

import argparse, math, random, sys
from collections import defaultdict

# --- sécurité affichage entiers géants ---
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

LOG10 = math.log(10.0)

# --------- utilitaires premiers / primalité ---------
_SMALL_PRIMES = [
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,
    211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331
]

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False
    # Miller–Rabin probabiliste (bases aléatoires)
    # suffisant pour de l'exploration
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    # 8 bases aléatoires + bases fixes pour fiabilité
    bases = [2, 3, 5, 7, 11, 13, 17]
    for _ in range(5):
        a = random.randrange(2, n - 2)
        bases.append(a)
    for a in bases:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        good = False
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                good = True
                break
        if not good:
            return False
    return True

def next_prime(n: int) -> int:
    if n <= 2:
        return 2
    p = n + 1 if n % 2 == 0 else n
    while True:
        if is_probable_prime(p):
            return p
        p += 2

# --------- Pollard–Rho (factorisation) ----------
def _rho_f(x, c, mod):
    return (x * x + c) % mod

def pollard_rho(n: int) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        x = random.randrange(2, n - 1)
        y = x
        c = random.randrange(1, n - 1)
        d = 1
        while d == 1:
            x = _rho_f(x, c, n)
            y = _rho_f(_rho_f(y, c, n), c, n)
            d = math.gcd(abs(x - y), n)
            if d == n:
                break
        if d > 1 and d < n:
            return d

def factorint(n: int, out: dict = None) -> dict:
    """Renvoie dict prime->exposant (factorisation complète)."""
    if out is None:
        out = defaultdict(int)
    if n < 2:
        return out
    for p in _SMALL_PRIMES:
        while n % p == 0:
            out[p] += 1
            n //= p
        if n == 1:
            return out
    if n == 1:
        return out
    if is_probable_prime(n):
        out[n] += 1
        return out
    d = pollard_rho(n)
    factorint(d, out)
    factorint(n // d, out)
    return out

# --------- σ(p^e) ----------
def sigma_p_pow(p: int, e: int) -> int:
    # σ(p^e) = (p^(e+1)-1)/(p-1)  exact
    return (pow(p, e + 1) - 1) // (p - 1)

# --------- outils D_structure ----------
def digits_lower_bound(primes_set) -> int:
    # borne dure : ceil(sum log10(p))
    s = 0.0
    for p in primes_set:
        s += math.log10(p)
    return int(math.ceil(s))

def digits_exact(primes_set, max_primes_for_exact=5000) -> int:
    # exact via produit, limité pour éviter les produits monstrueux
    if len(primes_set) > max_primes_for_exact:
        return None
    prod = 1
    for p in primes_set:
        prod *= p
    return len(str(prod))

# --------- génération d'ensembles initiaux ----------
def initial_scenario(y_start: int, scenario: str):
    y0 = y_start if is_probable_prime(y_start) else next_prime(y_start)
    if scenario.upper() == "A":
        # minimaliste : {p_min}
        return {y0}, None  # no special q
    elif scenario.upper() == "B":
        # "forme d'Euler" approximative : {q ≡ 1 mod 4, p1, p2}
        q = y0
        while q % 4 != 1:
            q = next_prime(q + 1)
        p1 = y0
        if p1 == q:
            p1 = next_prime(p1 + 1)
        p2 = next_prime(p1 + 1)
        return {q, p1, p2}, q
    else:
        raise ValueError("scenario doit être 'A' ou 'B'")

# --------- stratégie d'exposants ----------
def choose_exponent(p: int, special_q: int, strategy: str) -> int:
    st = strategy.lower()
    if st == "square":
        return 1 if (special_q is not None and p == special_q) else 2
    if st == "p_plus_1":
        return 1
    if st == "mixed":
        return 1 if (special_q is not None and p == special_q) else 2
    # défaut
    return 2

# --------- cœur de l’exploration ----------
def explore(y_start: int, scenario: str, max_depth: int, exp_strategy: str,
            respect_y_rude: bool, max_omega: int, exact_digits: bool, verbose: bool):
    known = set()
    frontier = set()
    init_set, special_q = initial_scenario(y_start, scenario)
    known |= init_set
    frontier |= init_set

    print("\n=== explore_npi_graph ===")
    print(f"y_start = {y_start}  | scenario = {scenario}  | exp_strategy = {exp_strategy}  | max_depth = {max_depth}")
    if special_q:
        print(f"special_q (≡1 mod 4) = {special_q}")
    print(f"respect_y_rude = {respect_y_rude}  | max_omega = {max_omega}  | exact_digits = {exact_digits}\n")

    def d_struct(prs):
        if exact_digits:
            de = digits_exact(prs)
            if de is not None:
                return de, True
        return digits_lower_bound(prs), False

    # profondeur 0 (état initial)
    d0, is_exact0 = d_struct(known)
    ratio0 = (d0 / max(1.0, len(known) * math.log10(max(2, y_start))))
    print("depth |  omega  |  D_structure  |  exact?  |  newest_prime  |  ratio / (ω log10 y)")
    newest = max(known)
    print(f"{0:5d} | {len(known):6d} | {d0:12d} | {str(is_exact0):7s} | {newest:13d} | {ratio0:16.6f}")

    for depth in range(1, max_depth + 1):
        if not frontier:
            if verbose:
                print(f"[stop] frontier vide à depth={depth-1}")
            break
        if len(known) >= max_omega:
            if verbose:
                print(f"[stop] ω atteint max_omega={max_omega}")
            break

        to_process = list(frontier)
        frontier = set()
        added_this_round = set()
        small_viols = 0

        for p in to_process:
            e = choose_exponent(p, special_q, exp_strategy)
            sig = sigma_p_pow(p, e)
            fac = factorint(sig)
            # ne garder que les clés (premiers)
            for r in fac.keys():
                if respect_y_rude and r < y_start:
                    small_viols += 1
                    continue
                if r not in known:
                    known.add(r)
                    added_this_round.add(r)
                    frontier.add(r)
                    if len(known) >= max_omega:
                        break
            if len(known) >= max_omega:
                break

        dval, is_exact = d_struct(known)
        newest = max(added_this_round) if added_this_round else max(known)
        ratio = (dval / max(1.0, len(known) * math.log10(max(2, y_start))))
        print(f"{depth:5d} | {len(known):6d} | {dval:12d} | {str(is_exact):7s} | {newest:13d} | {ratio:16.6f}")
        if verbose:
            print(f"   + nouveaux = {len(added_this_round)}  | violations y-rude (<y_start) ignorées = {small_viols}")

        if not added_this_round:
            if verbose:
                print(f"[stop] aucune expansion à depth={depth}")
            break

    print("\n[done]")
    return

# --------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Étape 5a : exploration de l'expansion forcée via σ(p^e)")
    ap.add_argument("--y-start", type=int, required=True, help="seuil p_min (y-rude)")
    ap.add_argument("--scenario", choices=["A","B"], default="B",
                    help="A: {p_min}; B: {q≡1 (mod 4), p1, p2}")
    ap.add_argument("--max-depth", type=int, default=6, help="nombre d'itérations (vagues)")
    ap.add_argument("--exp-strategy", choices=["square","p_plus_1","mixed"], default="square",
                    help="exposants e pour σ(p^e)")
    ap.add_argument("--respect-y-rude", action="store_true",
                    help="ignorer les facteurs r < y_start (sinon on les garde)")
    ap.add_argument("--max-omega", type=int, default=5000, help="garde-fou sur |F|")
    ap.add_argument("--exact-digits", action="store_true",
                    help="essaie un calcul exact des digits (limité, sinon borne dure)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    random.seed(42)
    explore(
        y_start=args.y_start,
        scenario=args.scenario,
        max_depth=args.max_depth,
        exp_strategy=args.exp_strategy,
        respect_y_rude=args.respect_y_rude,
        max_omega=args.max_omega,
        exact_digits=args.exact_digits,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()
