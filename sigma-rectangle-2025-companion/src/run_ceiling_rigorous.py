#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR

# =====================
# Réglages de précision
# =====================
getcontext().prec = 120  # haute précision pour logs géants
LN2 = Decimal(2).ln()
LN10 = Decimal(10).ln()

def dln(x: Decimal) -> Decimal:
    return x.ln()

def parse_big10power(s: str) -> int:
    s = s.strip().replace('_','')
    if '^' in s:
        base,exp = s.split('^',1)
        base = base.strip()
        exp = exp.strip()
        if base != '10':
            raise ValueError("Seul le format 10^D est supporté pour --x.")
        D = int(exp)
        return D + 1
    if not s.isdigit():
        raise ValueError("Pour --x, donnez '10^D' ou un entier décimal.")
    return len(s)

def compute_m_from_digits_y(D: int, y: int) -> (Decimal, Decimal, Decimal):
    y_dec = Decimal(y)
    ln_y = dln(y_dec)
    D_ln10_over_ln_y = (Decimal(D) * LN10) / ln_y
    # BUG FIX: utiliser la constante module-level ROUND_FLOOR
    m_dec = D_ln10_over_ln_y.to_integral_value(rounding=ROUND_FLOOR)
    return m_dec, ln_y, D_ln10_over_ln_y

def constant_ub_test(y: int, m: Decimal):
    y_dec = Decimal(y)
    z = Decimal(1) / (y_dec - 1)
    lhs_no = m * z
    if lhs_no < LN2:
        return ("NO", f"m*z < ln2 (lhs={lhs_no})")
    z2 = z*z
    lhs_ge = m * (z - z2/Decimal(2))
    if lhs_ge > LN2:
        return ("UNKNOWN_GE", f"m*(z - z^2/2) > ln2 (lhs={lhs_ge})")
    return ("UNKNOWN", "constant UB indéterminé (entre les deux tests)")

def solve_lnU_for_budget(y: int, ln_y: Decimal, m: Decimal, a: Decimal, b: Decimal, pad_factor: Decimal = Decimal('1.000001')) -> Decimal:
    y_dec = Decimal(y)
    C = y_dec / (ln_y - b)
    A = m + C
    if A <= 0:
        A = m if m > 1 else Decimal(2)
    lnA = dln(A)
    t = lnA + dln(lnA + a)
    for _ in range(50):
        t_next = lnA + dln(t + a)
        if abs(t_next - t) < Decimal('1e-30'):
            t = t_next
            break
        t = t_next
    return t + dln(pad_factor)

def sandwich_decide(y: int, ln_y: Decimal, m: Decimal, a_pi: Decimal, b_pi: Decimal):
    t = solve_lnU_for_budget(y, ln_y, m, a_pi, b_pi)
    eps = Decimal('1e-18')
    denom_up = (ln_y - 1) if (ln_y - 1) > eps else (ln_y + eps)
    denom_lo = (ln_y + 1)
    upper_ratio = (t + 1) / denom_up
    lower_ratio = (t - 1) / denom_lo
    two = Decimal(2)
    if upper_ratio < two:
        return ("likely_NO", t, upper_ratio, lower_ratio)
    if lower_ratio > two:
        return ("likely_YES", t, upper_ratio, lower_ratio)
    return ("UNKNOWN", t, upper_ratio, lower_ratio)

def main():
    ap = argparse.ArgumentParser(description="Abundance ceiling (rigorous + sandwich)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--digits", type=int, help="D = digits(x)")
    g.add_argument("--x", type=str, help="x as '10^D' or big integer literal")

    ap.add_argument("--y", type=int, required=True, help="minimal prime factor constraint (integer ≥ 3)")

    ap.add_argument("--mertens-upper", action="store_true", help="(NO-only) hook conservé")
    ap.add_argument("--mertens-lower", action="store_true", help="(YES-only) hook conservé")

    ap.add_argument("--sandwich", action="store_true", help="fast sandwich (no primes), likely YES/NO if clear")
    ap.add_argument("--pi-a", type=float, default=1.0, help="π(x) ≥ x/(ln x + a), a≥0 (default 1.0)")
    ap.add_argument("--pi-b", type=float, default=1.0, help="π(x) ≤ x/(ln x − b), b≥0 (default 1.0)")

    args = ap.parse_args()

    if args.digits is not None:
        D = int(args.digits)
    else:
        D = parse_big10power(args.x)

    y = int(args.y)
    if y < 3:
        raise SystemExit("y doit être ≥ 3")

    m, ln_y, _ = compute_m_from_digits_y(D, y)

    print("=== Abundance ceiling (rigorous + sandwich) ===")
    print(f"digits(x) = {D}")
    print(f"y = {y}")
    print(f"m = floor(log_y x) = {m}")

    tag, info = constant_ub_test(y, m)
    print("\n[Constant UB]  (rigorous)")
    if tag == "NO":
        print(f"  decision: NO (provable)")
        print(f"  reason  : {info}")
        if not args.sandwich:
            return
    elif tag == "UNKNOWN_GE":
        print(f"  decision: UNKNOWN")
        print(f"  note    : {info} ⇒ (y/(y-1))^m ≥ 2, mais ce test ne prouve pas YES.")
    else:
        print(f"  decision: UNKNOWN")
        print(f"  note    : {info}")

    if args.sandwich:
        with localcontext() as ctx:
            ctx.prec = 120
            a_pi = Decimal(str(args.pi_a))
            b_pi = Decimal(str(args.pi_b))
            s_tag, t_lnU, up, lo = sandwich_decide(y, ln_y, m, a_pi, b_pi)

        print("\n[Sandwich Mertens-like]  (fast, no primes)")
        print(f"  ln U (solved) ≈ {t_lnU}")
        print(f"  upper_ratio ≈ {up}")
        print(f"  lower_ratio ≈ {lo}")
        if s_tag == "likely_NO":
            print("  decision: likely NO")
            print("  note    : upper_ratio < 2 (conservateur).")
        elif s_tag == "likely_YES":
            print("  decision: likely YES")
            print("  note    : lower_ratio > 2 (conservateur).")
        else:
            print("  decision: UNKNOWN")
            print("  note    : sandwich trop serré; envisage un stream partiel ou des bornes π plus fines.")

    print("\nNote:")
    print("  - 'NO (provable)' vient du test constant UB (m*z < ln2).")
    print("  - Le sandwich rend 'likely YES/NO' sans énumérer de premiers; rapide mais pas une preuve universelle.")
    print("  - Le GPU n’aide pas ici: on fait de l’arithmétique haute précision sérielle (Decimal).")
    print("  - Pour un YES rigoureux sans primes, il faut des bornes Mertens explicites plus serrées ou un stream partiel ∑1/p avec erreur contrôlée.")
if __name__ == "__main__":
    main()
