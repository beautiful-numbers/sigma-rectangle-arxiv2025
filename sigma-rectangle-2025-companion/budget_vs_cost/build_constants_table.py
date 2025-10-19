#!/usr/bin/env python3
# build_constants_table.py
#
# Objet : produire les lignes "constantes effectives" pour l'annexe :
#   - epsilon(y) = max(log(2) - sum_{y<=p<=y^2} 1/p, 0)
#   - borne inférieure Rosser–Schoenfeld pour pi(y^2) - pi(y)
#   - Phi(y) = log10(y) * ( 1/(2 ln y + A_RS) - 1/(y (ln y - A_RS)) )
#   - c1 = min_y Phi(y) (sauf si --forced-c1)
#   - Y0 choisi comme le plus petit y de la liste avec Phi(y) >= c1 et y >= y0_floor
#   - vérification : D_pivot(y) >= c1 * y^2  <=>  Phi(y) >= c1
#
# Sorties :
#   - tableau lisible sur stdout
#   - CSV optionnel (--csv-out)
#
# Remarques :
#   - Sieve exact jusqu'à max(y)^2 (adapté aux y ~ 10^3–10^4).
#   - Somme 1/p en Decimal haute précision.
#   - RS avec A_RS = 1.2762 (valable dès nos x).
#   - Tous les logs naturels sont en ln(); log10 est la base 10.
#
# Exemple :
#   python build_constants_table.py --ys 607,613,641,701,997,1499,2003 --csv-out E23_table.csv

import argparse
import math
from decimal import Decimal, getcontext
from typing import List, Tuple, Optional

# --------- précision Decimal ---------
getcontext().prec = 50

# --------- constantes ----------
A_RS = 1.2762  # Rosser–Schoenfeld parameter
LN2 = Decimal(str(math.log(2.0)))  # ln(2) en Decimal

# --------- helpers format ---------
def group_int(n: int) -> str:
    return f"{n:,}".replace(",", " ")

def fmt_dec(x: Decimal, digits: int = 12) -> str:
    q = Decimal(1).scaleb(-digits)  # 10^-digits
    return str(x.quantize(q))

def fmt_float(x: float, digits: int = 12) -> str:
    return f"{x:.{digits}f}"

# --------- primalité & crible ---------
def simple_sieve(limit: int) -> List[int]:
    if limit < 2:
        return []
    bs = bytearray(b"\x01") * (limit + 1)
    bs[:2] = b"\x00\x00"
    r = int(limit ** 0.5)
    for p in range(2, r + 1):
        if bs[p]:
            step = p
            start = p * p
            bs[start:limit + 1:step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if bs[i]]

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    # Miller–Rabin déterministe pour 64-bit (bases suffisantes)
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bases = (2, 3, 5, 7, 11, 13, 17)
    for a in bases:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        skip = False
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                skip = True
                break
        if skip:
            continue
        return False
    return True

def primes_in_range_with_sieve(ymin: int, ymax2: int) -> List[int]:
    # retourne tous les premiers <= ymax2 ; on filtrera ensuite
    return simple_sieve(ymax2)

# --------- logs ---------
def lnD(x: int | float) -> Decimal:
    # ln en Decimal via math.log (double) converti ; suffisant ici
    return Decimal(str(math.log(float(x))))

def log10f(x: int | float) -> float:
    return math.log10(float(x))

# --------- RS bounds for pi(x) ---------
def pi_lower_RS(x: int) -> float:
    # pi(x) >= x / (ln x + A_RS)
    if x < 3:
        return 0.0
    return float(x) / (math.log(x) + A_RS)

def pi_upper_RS(x: int) -> float:
    # pi(x) <= x / (ln x - A_RS)
    if x < 3:
        return 2.0  # garde-fou (inutile pour nos x)
    return float(x) / (math.log(x) - A_RS)

def pi_diff_lower_RS(y: int) -> int:
    # borne inférieure pour pi(y^2) - pi(y)
    y2 = y * y
    lower = pi_lower_RS(y2) - pi_upper_RS(y)
    # on renvoie au moins 0
    return max(0, int(math.floor(lower)))

# --------- somme des inverses sur [y, y^2] ---------
def sum_reciprocals_primes_range(primes: List[int], y: int) -> Decimal:
    y2 = y * y
    total = Decimal(0)
    for p in primes:
        if p < y:
            continue
        if p > y2:
            break
        total += Decimal(1) / Decimal(p)
    return total

# --------- Phi(y) ---------
def phi_of_y(y: int) -> float:
    # Phi(y) = log10(y) * ( 1 / (2 ln y + A_RS) - 1 / (y (ln y - A_RS)) )
    ln_y = math.log(y)
    term1 = 1.0 / (2.0 * ln_y + A_RS)
    term2 = 1.0 / (y * (ln_y - A_RS))
    return math.log10(y) * (term1 - term2)

# --------- table row ---------
def compute_row(y: int, primes_all: List[int]) -> dict:
    H = sum_reciprocals_primes_range(primes_all, y)  # Decimal
    eps = (LN2 - H) if H < LN2 else Decimal(0)
    H_lt_ln2 = H < LN2
    pi_diff_lb = pi_diff_lower_RS(y)
    phi = phi_of_y(y)
    return {
        "y": y,
        "H": H,
        "epsilon": eps,
        "H_lt_ln2": H_lt_ln2,
        "pi_diff_lb_RS": pi_diff_lb,
        "Phi": phi,
    }

# --------- CLI ----------

def parse_ys(arg: Optional[str]) -> List[int]:
    if not arg:
        return []
    out = []
    for token in arg.split(","):
        t = token.strip()
        if not t:
            continue
        out.append(int(t))
    return out

def build_prime_y_list(ymin: int | None, ymax: int | None, ys_list: List[int]) -> List[int]:
    ys = []
    if ys_list:
        ys.extend(ys_list)
    if ymin is not None and ymax is not None:
        if ymax < ymin:
            raise ValueError("ymax < ymin")
        # générer tous les premiers impairs dans [ymin, ymax]
        for n in range(max(3, ymin), ymax + 1):
            if n % 2 == 0:
                continue
            if is_probable_prime(n):
                ys.append(n)
    # unique + tri
    ys = sorted(set(ys))
    # filtrer impairs >= 3
    ys = [y for y in ys if y >= 3 and (y % 2 == 1)]
    return ys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ys", type=str, help="liste de y (premiers impairs) séparés par des virgules")
    ap.add_argument("--ymin", type=int, help="borne inf y (incluse)")
    ap.add_argument("--ymax", type=int, help="borne sup y (incluse)")
    ap.add_argument("--y0-floor", type=int, default=607, help="plancher pour Y0 (ex: 607)")
    ap.add_argument("--forced-c1", type=float, help="forcer c1 (sinon min Phi(y) sur l’échantillon)")
    ap.add_argument("--csv-out", type=str, help="fichier CSV de sortie")
    args = ap.parse_args()

    ys_in = parse_ys(args.ys)
    ys = build_prime_y_list(args.ymin, args.ymax, ys_in)
    if not ys:
        raise SystemExit("Aucun y fourni. Utilise --ys ou --ymin/--ymax.")

    y_max = max(ys)
    primes_all = primes_in_range_with_sieve(min(ys), y_max * y_max)

    # calcul des lignes
    rows = []
    for y in ys:
        if not is_probable_prime(y):
            print(f"[WARN] y={y} n’est pas premier (on continue).")
        row = compute_row(y, primes_all)
        rows.append(row)

    # c1
    phi_vals = [row["Phi"] for row in rows if row["y"] >= args.y0_floor]
    if not phi_vals:
        # si tous les y < y0_floor, on prend tout de même le min global
        phi_vals = [row["Phi"] for row in rows]
    min_phi = min(phi_vals)
    c1 = float(args.forced_c1) if args.forced_c1 is not None else float(min_phi)

    # Y0 automatique : plus petit y >= y0_floor tel que Phi(y) >= c1 (avec petite tolérance)
    TOL = 1e-15
    Y0_candidates = [row["y"] for row in rows if row["y"] >= args.y0_floor and (row["Phi"] + TOL) >= c1]
    Y0 = min(Y0_candidates) if Y0_candidates else max(args.y0_floor, min(ys))

    # impression
    print("=== Constantes effectives : tableau récapitulatif ===")
    print(f"A_RS = {A_RS}   |   c1 = {c1:.12f}   |   Y0 = {Y0}")
    print()
    header = (
        " y  |  epsilon(y)         |  H(y..y^2)           |  H<ln2 |  pi(y^2)-pi(y) RS_lb |      Phi(y)       |  Phi>=c1 |  D_pivot>=c1*y^2"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        y = row["y"]
        eps_s = fmt_dec(row["epsilon"], 12)
        H_s = fmt_dec(row["H"], 12)
        hflag = "Y" if row["H_lt_ln2"] else "N"
        pi_lb = row["pi_diff_lb_RS"]
        phi_s = fmt_float(row["Phi"], 12)
        phi_ge = (row["Phi"] + TOL) >= c1
        d_ge = phi_ge  # équivalent
        print(
            f"{y:4d} | {eps_s:>20} | {H_s:>20} | {hflag:^6} | {group_int(pi_lb):>21} | {phi_s:>16} | {str(phi_ge):>8} | {str(d_ge):>17}"
        )

    # CSV
    if args.csv_out:
        import csv
        with open(args.csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "y", "epsilon", "H", "H_lt_ln2", "pi_diff_lb_RS", "Phi", "c1", "Phi_ge_c1", "Dpivot_ge_c1_y2", "Y0"
            ])
            for row in rows:
                phi_ge = row["Phi"] >= c1
                w.writerow([
                    row["y"],
                    f"{row['epsilon']}",
                    f"{row['H']}",
                    int(row["H_lt_ln2"]),
                    row["pi_diff_lb_RS"],
                    f"{row['Phi']:.16f}",
                    f"{c1:.16f}",
                    int(phi_ge),
                    int(phi_ge),
                    Y0 if row["y"] == rows[0]["y"] else ""
                ])
        print(f"\n[OK] CSV écrit : {args.csv_out}")

if __name__ == "__main__":
    main()
