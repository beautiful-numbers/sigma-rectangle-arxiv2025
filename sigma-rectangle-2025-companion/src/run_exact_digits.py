#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import sys

# lever toute limite de conversion d'entiers géants en str
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

def ndigits(n: int) -> int:
    return len(str(n)) if n != 0 else 1

def digits_of_power(y: int, m: int) -> int:
    return ndigits(pow(y, m))

def find_m_with_exact_digits(y: int, D: int):
    if y < 2: 
        return None
    m0 = max(1, int(math.ceil((D - 1) / math.log10(y))))
    # petite fenêtre autour de m0
    for delta in range(-50, 51):
        m = m0 + delta
        if m > 0 and digits_of_power(y, m) == D:
            return m
    return None

def first_m_with_at_least_D_digits(y: int, D: int) -> int:
    m = max(1, int(math.floor((D - 1) / math.log10(y))))
    while digits_of_power(y, m) < D:
        m += 1
    while m > 1 and digits_of_power(y, m - 1) >= D:
        m -= 1
    return m

def main():
    ap = argparse.ArgumentParser(description="Génère l’ancre y^m à D chiffres (si elle existe) et la frontière budget (y^m, y^m-1).")
    ap.add_argument("--y", type=int, required=True, help="Base (ex: un premier).")
    ap.add_argument("--digits", type=int, required=True, help="Nombre de chiffres D visé pour y^m.")
    ap.add_argument("--print-numbers", dest="print_numbers", action="store_true",
                    help="Imprimer aussi les grands entiers.")
    args = ap.parse_args()

    y = args.y
    D = args.digits

    print(f"\n=== run_exact_digits ===")
    print(f"y = {y}")
    print(f"D = {D} (nombre de chiffres visé pour y^m)")

    m_exact = find_m_with_exact_digits(y, D)
    if m_exact is not None:
        N_hi = pow(y, m_exact)
        N_lo = N_hi - 1
        d_hi = ndigits(N_hi)
        d_lo = ndigits(N_lo)
        print("\n[Exact match trouvé]")
        print(f"  m = {m_exact}  avec  digits(y^m) = {d_hi} (== D)")
        print(f"  N_hi = y^m       (palier budget)")
        print(f"  N_lo = y^m - 1   (voisin immédiat en dessous)")
        print(f"  digits(N_hi) = {d_hi},  digits(N_lo) = {d_lo}")
        if args.print_numbers:
            print("\nN_hi =")
            print(N_hi)
            print("\nN_lo =")
            print(N_lo)
    else:
        print("\n[Aucun m ne donne exactement D chiffres pour y^m]")
        m_up = first_m_with_at_least_D_digits(y, D)
        N_border = pow(y, m_up)
        N_border_minus = N_border - 1
        d_b = ndigits(N_border)
        d_bm = ndigits(N_border_minus)
        print(f"  Frontière budget (palier D) : m_up = {m_up}")
        print(f"  y^m_up a {d_b} chiffres (>= D), (y^m_up - 1) en a {d_bm}.")
        if args.print_numbers:
            print("\n(y^m_up) =")
            print(N_border)
            print("\n(y^m_up - 1) =")
            print(N_border_minus)

    # Toujours proposer l’ancre “exact D chiffres” universelle
    TD = pow(10, D) - 1
    print(f"\n[Ancre exacte à D chiffres]")
    print(f"  T_D = 10^{D} - 1  (a exactement {D} chiffres)")
    if args.print_numbers:
        print("\nT_D =")
        print(TD)

if __name__ == "__main__":
    main()
