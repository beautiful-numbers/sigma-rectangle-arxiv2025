#!/usr/bin/env python3
# calculate_d_structure.py
# D_structure(y, ω) = digits( ∏_{j=1..ω} p_j ) where p_j are the ω smallest primes ≥ y
# - Génération des premiers : crible si y ≤ 1e7, sinon Miller–Rabin + next_prime
# - Produit exact en grands entiers
# - Nombre de chiffres exact via recherche binaire sur 10^k (pas de str géant)

import argparse, math, sys

# Lever toute limite de conversion si jamais on imprime un entier (nous n'en avons pas besoin ici,
# mais on met la protection au cas où tu ajoutes un --print-product plus tard)
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

LOG10 = math.log(10.0)

# ---------- outils premiers ----------
def simple_sieve(limit: int):
    bs = bytearray(b"\x01") * (limit + 1)
    bs[:2] = b"\x00\x00"
    r = int(limit**0.5)
    for p in range(2, r + 1):
        if bs[p]:
            step = p
            start = p * p
            bs[start:limit+1:step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if bs[i]]

_MR_BASES_64 = (2, 3, 5, 7, 11, 13, 17)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2,3,5,7,11,13,17,19,23,29)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    # écrire n-1 = d * 2^s
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    for a in _MR_BASES_64:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
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

def primes_from_y(y: int, m: int):
    """
    Génère les m plus petits nombres premiers ≥ y.
    - Crible si y ≤ 1e7 avec une borne raisonnable
    - Sinon Miller–Rabin/next_prime
    """
    if m <= 0:
        return []
    if y <= 10_000_000:
        ly = math.log(max(3, y))
        bound = int(y + m * (ly + math.log(ly + 1)))
        bound = max(bound, y + 2000)
        primes = simple_sieve(bound)
        # couper < y
        import bisect
        i = bisect.bisect_left(primes, y)
        out = primes[i:i+m]
        if len(out) < m:
            # compléter en MR au besoin
            p = primes[-1] + 1 if primes else y
            if p % 2 == 0:
                p += 1
            while len(out) < m:
                if is_probable_prime(p):
                    out.append(p)
                p += 2
        return out
    else:
        out = []
        p = y if is_probable_prime(y) else next_prime(y)
        out.append(p)
        while len(out) < m:
            p = next_prime(p + 1)
            out.append(p)
        return out

# ---------- digits exacts par recherche binaire ----------
def digits_bigint_exact(n: int, lower_hint: int = 1, upper_hint: int | None = None) -> int:
    """
    Renvoie digits(n) = 1 + floor(log10 n) sans conversion en chaîne.
    On cherche le plus petit k tel que 10^k > n.

    lower_hint : borne inférieure plausible (ex: ceil(ω log10 y))
    upper_hint : borne supérieure (ex: int(sum log10 p) + marge)
    """
    if n == 0:
        return 1
    if n < 10:
        return 1

    # s'assurer d'une borne basse ≥ 1
    low = max(1, int(lower_hint))
    # trouver une borne haute
    if upper_hint is None:
        high = low
        ten_pow = pow(10, high)
        while ten_pow <= n:
            high *= 2
            ten_pow = pow(10, high)
    else:
        high = max(low + 1, int(upper_hint))
        ten_pow = pow(10, high)
        while ten_pow <= n:
            high *= 2
            ten_pow = pow(10, high)

    # maintenant on a 10^high > n (ten_pow défini)
    # garantir aussi 10^low ≤ n
    while pow(10, low) > n:
        low = max(1, low - 1)

    # binaire classique sur k
    # invariant: 10^low ≤ n < 10^high
    while low + 1 < high:
        mid = (low + high) // 2
        if pow(10, mid) <= n:
            low = mid
        else:
            high = mid
    return high

# ---------- calcul principal ----------
def main():
    ap = argparse.ArgumentParser(description="Compute D_structure(y, ω) with exact big-int product of ω smallest primes ≥ y")
    ap.add_argument("--y", type=int, required=True, help="plus petit premier autorisé (p_min ≥ y)")
    ap.add_argument("--omega", type=int, required=True, help="nombre de facteurs premiers distincts ω")
    ap.add_argument("--print-primes", action="store_true", help="afficher la liste des ω premiers utilisés (attention: long)")
    args = ap.parse_args()

    y = int(args.y)
    w = int(args.omega)
    if y < 2 or w < 0:
        print("Paramètres invalides.")
        return
    if y == 2:
        y = 3  # on veut des impairs pour coller au cadre NPI, mais l'outil fonctionne pour tout y

    print("\n=== calculate_d_structure ===")
    print(f"y = {y}")
    print(f"omega = {w}")

    # borne inférieure immédiate (produit ≥ y^ω)
    digits_lb = int(math.ceil(w * math.log10(y))) if w > 0 else 1

    # générer les ω premiers ≥ y
    plist = primes_from_y(y, w)
    if args.print_primes:
        print(f"\nprimes (ω={w}) >= {y} :")
        print(plist)

    if w == 0:
        print("\n--- Résultat ---")
        print("Produit vide = 1")
        print("D_structure = 1")
        print(f"lower_bound = {digits_lb}")
        return

    last_p = plist[-1]
    # estimation via somme des log10 (utile pour borner la recherche binaire)
    sum_log10 = 0.0
    for p in plist:
        sum_log10 += math.log10(p)
    digits_est = int(sum_log10) + 1  # estimation haute (souvent exacte), pas certifiée

    # produit exact
    P = 1
    for p in plist:
        P *= p

    # digits exacts par binaire, bornes: [digits_lb, digits_est + marge]
    digits_exact = digits_bigint_exact(P, lower_hint=digits_lb, upper_hint=digits_est + 5)

    print("\n--- Résultat ---")
    print(f"last_prime_used = {last_p}")
    print(f"D_structure (exact) = {digits_exact}")
    print(f"lower_bound = ceil(ω * log10(y)) = {digits_lb}")
    print(f"est_digits_from_logs ≈ {digits_est}  (non-certifié, utilisé comme borne)")
    print(f"check: exact ≥ lower_bound ? {'YES' if digits_exact >= digits_lb else 'NO'}")

if __name__ == "__main__":
    main()
