#!/usr/bin/env python3
# calculate_prime_density.py
# Affiche, tranche par tranche, le nombre de nombres premiers dans [start, end]
# en utilisant un crible segmenté. Pas de CSV : sortie directe dans la console.

import argparse, math

def simple_sieve(limit: int):
    if limit < 2:
        return []
    bs = bytearray(b"\x01") * (limit + 1)
    bs[:2] = b"\x00\x00"
    r = int(limit**0.5)
    for p in range(2, r + 1):
        if bs[p]:
            step = p
            start = p * p
            bs[start:limit+1:step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if bs[i]]

def count_primes_segment(a: int, b: int, base_primes):
    """Compte les premiers dans [a, b] via crible segmenté."""
    n = b - a + 1
    if n <= 0:
        return 0
    isprime = bytearray(b"\x01") * n

    # 0 et 1 ne sont pas premiers
    if a <= 1 <= b:
        isprime[1 - a] = 0
    if a == 0:
        isprime[0] = 0

    for p in base_primes:
        pp = p * p
        if pp > b:
            break
        start = pp if pp >= a else ((a + p - 1) // p) * p
        # ne pas marquer p lui-même si p ∈ [a, b]
        if start == p:
            start += p
        for m in range(start, b + 1, p):
            isprime[m - a] = 0

    return int(sum(isprime))

def approx_interval_count(a: int, b: int):
    """Approximation PNT simple: π(b)-π(a-1) ≈ b/ln b - (a-1)/ln(a-1)."""
    def term(x):
        if x < 3:
            return 0.0
        return x / math.log(x)
    return term(b) - term(a - 1)

def run(start: int, end: int, slice_size: int, show_theory: bool, show_density: bool):
    if start > end:
        start, end = end, start
    if slice_size <= 0:
        raise ValueError("slice-size doit être > 0")

    # Prépare les premiers jusqu'à sqrt(end) pour le crible segmenté
    base_limit = int(math.isqrt(end)) + 1
    base_primes = simple_sieve(base_limit)

    total_count = 0
    total_len = 0
    idx = 0

    print("\n=== Prime density by slice ===")
    print(f"interval = [{start}, {end}]  |  slice_size = {slice_size}\n")

    for a in range(start, end + 1, slice_size):
        b = min(a + slice_size - 1, end)
        idx += 1
        cnt = count_primes_segment(a, b, base_primes)
        length = b - a + 1
        total_count += cnt
        total_len += length

        if show_theory:
            approx = approx_interval_count(a, b)
            err = cnt - approx
            rel = (err / cnt) if cnt else float('inf')
            if show_density:
                dens = cnt / length
                print(f"slice {idx:2d}: [{a:,}..{b:,}]  len={length:,}  primes={cnt:,}  "
                      f"density={dens:.6f}  approx={approx:,.2f}  err={err:,.2f}  rel={rel:.3%}")
            else:
                print(f"slice {idx:2d}: [{a:,}..{b:,}]  len={length:,}  primes={cnt:,}  "
                      f"approx={approx:,.2f}  err={err:,.2f}  rel={rel:.3%}")
        else:
            if show_density:
                dens = cnt / length
                print(f"slice {idx:2d}: [{a:,}..{b:,}]  len={length:,}  primes={cnt:,}  density={dens:.6f}")
            else:
                print(f"slice {idx:2d}: [{a:,}..{b:,}]  len={length:,}  primes={cnt:,}")

    print("\n--- Totaux ---")
    if show_density:
        print(f"total_len={total_len:,}  total_primes={total_count:,}  density={total_count/total_len:.6f}")
    else:
        print(f"total_len={total_len:,}  total_primes={total_count:,}")
    if show_theory:
        approx_tot = approx_interval_count(start, end)
        err_tot = total_count - approx_tot
        rel_tot = (err_tot / total_count) if total_count else float('inf')
        print(f"approx_total={approx_tot:,.2f}  err_total={err_tot:,.2f}  rel_total={rel_tot:.3%}")

def main():
    ap = argparse.ArgumentParser(description="Compteur de nombres premiers par tranches (crible segmenté).")
    ap.add_argument("--start", type=int, required=True, help="début de l'intervalle inclus")
    ap.add_argument("--end", type=int, required=True, help="fin de l'intervalle inclus")
    ap.add_argument("--slice-size", type=int, required=True, help="taille de chaque tranche")
    ap.add_argument("--theory", action="store_true", help="affiche aussi l'approximation PNT et l'erreur")
    ap.add_argument("--density", action="store_true", help="affiche la densité (primes/longueur) pour chaque tranche")
    args = ap.parse_args()
    run(args.start, args.end, args.slice_size, args.theory, args.density)

if __name__ == "__main__":
    main()
