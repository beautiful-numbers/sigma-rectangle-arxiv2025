#!/usr/bin/env python3
# calculate_m_requis.py
# Usage:
#   python calculate_m_requis.py --y 79                # résumé compact (recommandé)
#   python calculate_m_requis.py --y 79 --save-json cert_y79.json   # export complet JSON
#   python calculate_m_requis.py --y 79 --reduced --save-json cert_y79_reduced.json
#
# Sortie console = certificat strict mais court:
#   - y, m_requis, p_m, D_block = floor(m*log10(y))
#   - Vérifs EXACTES: UB_{m-1} < 2  et  UB_m > 2
#   - Empreintes SHA-256 et previews hex (num/den) de U_prev et U_curr
#   - Hash de la liste W_y (les m plus petits premiers >= y), et premiers de tête/fin
#
# Option --save-json : écrit un certificat complet (avec fractions, éventuellement réduites)
#                      pour archivage ou vérif externe.
# Aucune dépendance externe. Tout est en arithmétique entière exacte.

import argparse, math, sys, json, hashlib

# Lever d'éventuelles limites d'impression des très grands entiers (si jamais on exporte en JSON)
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

LOG10 = math.log(10.0)
LOG10_2 = math.log10(2.0)

# ---------- primalité (Miller–Rabin déterministe < 2^64) ----------
_MR_BASES_64 = (2, 3, 5, 7, 11, 13, 17)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = [2,3,5,7,11,13,17,19,23,29]
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

def prime_iter_from(y: int):
    p = y if is_probable_prime(y) else next_prime(y)
    while True:
        yield p
        p = next_prime(p + 1)

# ---------- utilitaires d’empreinte ----------
def int_sha256_hex(n: int) -> str:
    if n == 0:
        b = b"\x00"
    else:
        b = n.to_bytes((n.bit_length() + 7) // 8, "big", signed=False)
    return hashlib.sha256(b).hexdigest()

def hex_preview(n: int, k: int = 16) -> str:
    if n == 0:
        return "00"
    h = n.to_bytes((n.bit_length() + 7) // 8, "big").hex()
    return h if len(h) <= 2*k else (h[:k] + "..." + h[-k:])

def digits10_exact(n: int) -> int:
    # calcule le nombre de chiffres en base 10 sans convertir toute la valeur en chaîne
    if n == 0:
        return 1
    d = int((n.bit_length() - 1) * 0.30102999566398114) + 1  # ~ log10(2)
    # correction locale
    p = pow(10, d - 1)
    if n < p:
        while n < p:
            d -= 1
            p //= 10
    else:
        while True:
            q = p * 10
            if n >= q:
                d += 1
                p = q
            else:
                break
    return d

# ---------- cœur : calcul de m_requis(y) et certificat ----------
def compute_m_requis_and_certificate(y: int, keep_primes: bool):
    assert y >= 3
    primes = [] if keep_primes else None

    num = 1   # ∏ p
    den = 1   # ∏ (p-1)
    m = 0
    it = prime_iter_from(y)
    p_last = None

    while True:
        p = next(it)
        p_last = p
        new_num = num * p
        new_den = den * (p - 1)

        # franchit-on 2 strictement ?
        # On veut strictement > 2, pas >=
        if new_num > 2 * new_den:
            # m_requis = m+1 ; U_prev = (num/den) ; U_curr = (new_num/new_den)
            m_requis = m + 1
            if keep_primes:
                primes.append(p)
            U_prev_num, U_prev_den = num, den
            U_curr_num, U_curr_den = new_num, new_den
            break

        # sinon, on empile et on continue
        num, den = new_num, new_den
        m += 1
        if keep_primes:
            primes.append(p)

    # vérifs strictes et exactes
    assert U_prev_num * 1 < 2 * U_prev_den, "UB_{m-1} doit être strictement < 2"
    assert U_curr_num * 1 > 2 * U_curr_den, "UB_m doit être strictement > 2"

    return {
        "y": y,
        "m": m_requis,
        "p_m": p_last,
        "W_primes": primes,  # None si non demandé
        "U_prev_num": U_prev_num,
        "U_prev_den": U_prev_den,
        "U_curr_num": U_curr_num,
        "U_curr_den": U_curr_den,
    }

def reduce_fraction(num: int, den: int):
    import math
    g = math.gcd(num, den)
    if g != 1:
        num //= g
        den //= g
    return num, den

# ---------- impression compacte ----------
def print_compact_summary(res: dict, reduced: bool):
    y = res["y"]; m = res["m"]; p_m = res["p_m"]
    upn, upd = res["U_prev_num"], res["U_prev_den"]
    ucn, ucd = res["U_curr_num"], res["U_curr_den"]

    if reduced:
        upn, upd = reduce_fraction(upn, upd)
        ucn, ucd = reduce_fraction(ucn, ucd)

    D_block = int(math.floor(m * math.log10(y)))

    # empreintes
    upn_sha = int_sha256_hex(upn); upd_sha = int_sha256_hex(upd)
    ucn_sha = int_sha256_hex(ucn); ucd_sha = int_sha256_hex(ucd)

    print("\n=== calculate_m_requis (compact certificate) ===")
    print(f"y = {y}")
    print(f"m_requis(y) = {m}")
    print(f"p_m (dernier premier empilé) = {p_m}")
    print(f"D_block = floor(m * log10(y)) = {D_block}")
    print("\nVérifs exactes (strictes):")
    print(f"  UB_{{m-1}}(y) < 2  : {'YES' if upn * 1 < 2 * upd else 'NO'}")
    print(f"  UB_m(y)   > 2  : {'YES' if ucn * 1 > 2 * ucd else 'NO'}")

    print("\nU_prev = UB_{m-1}(y)  (empreintes)")
    print(f"  num_digits = {digits10_exact(upn)}  den_digits = {digits10_exact(upd)}")
    print(f"  num_sha256 = {upn_sha}")
    print(f"  den_sha256 = {upd_sha}")
    print(f"  num_hex = {hex_preview(upn)}")
    print(f"  den_hex = {hex_preview(upd)}")

    print("\nU_curr = UB_{m}(y)  (empreintes)")
    print(f"  num_digits = {digits10_exact(ucn)}  den_digits = {digits10_exact(ucd)}")
    print(f"  num_sha256 = {ucn_sha}")
    print(f"  den_sha256 = {ucd_sha}")
    print(f"  num_hex = {hex_preview(ucn)}")
    print(f"  den_hex = {hex_preview(ucd)}")

    # hash de la liste de premiers (si disponible)
    if res["W_primes"] is not None:
        # on hashe la concaténation binaire des entiers pour éviter une très longue chaîne
        h = hashlib.sha256()
        for p in res["W_primes"]:
            b = p.to_bytes((p.bit_length()+7)//8, "big", signed=False)
            # longueur + contenu pour éviter collisions triviales
            h.update(len(b).to_bytes(4, "big"))
            h.update(b)
        w_hash = h.hexdigest()

        head = res["W_primes"][:8]
        tail = res["W_primes"][-8:] if len(res["W_primes"]) > 8 else []
        print("\nW_y (les m plus petits premiers ≥ y)")
        print(f"  count = {len(res['W_primes'])}, sha256 = {w_hash}")
        print(f"  head = {head}")
        if tail:
            print(f"  tail = {tail}")

# ---------- export JSON (optionnel) ----------
def save_full_json(path: str, res: dict, reduced: bool):
    upn, upd = res["U_prev_num"], res["U_prev_den"]
    ucn, ucd = res["U_curr_num"], res["U_curr_den"]
    if reduced:
        upn, upd = reduce_fraction(upn, upd)
        ucn, ucd = reduce_fraction(ucn, ucd)

    obj = {
        "y": res["y"],
        "m": res["m"],
        "p_m": res["p_m"],
        "W_primes": res["W_primes"],  # peut être None
        "U_prev": {
            "num": str(upn),
            "den": str(upd),
            "num_digits": digits10_exact(upn),
            "den_digits": digits10_exact(upd),
            "sha256_num": int_sha256_hex(upn),
            "sha256_den": int_sha256_hex(upd),
        },
        "U_curr": {
            "num": str(ucn),
            "den": str(ucd),
            "num_digits": digits10_exact(ucn),
            "den_digits": digits10_exact(ucd),
            "sha256_num": int_sha256_hex(ucn),
            "sha256_den": int_sha256_hex(ucd),
        },
        "checks": {
            "U_prev_lt_2": (upn * 1 < 2 * upd),
            "U_curr_gt_2": (ucn * 1 > 2 * ucd),
        },
        "D_block": int(math.floor(res["m"] * math.log10(res["y"]))),
        "reduced": bool(reduced),
        "version": "1.2",
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False)
    print(f"\n[JSON] certificat écrit dans: {path}")

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Certificat strict pour m_requis(y): UB_{m-1}(y) < 2 < UB_m(y)")
    ap.add_argument("--y", type=int, required=True, help="borne minimale des premiers (p_min ≥ y)")
    ap.add_argument("--reduced", action="store_true", help="réduire U_prev/U_curr (gcd) dans l’export JSON et les empreintes")
    ap.add_argument("--save-json", type=str, help="chemin de sortie JSON (export complet)")
    ap.add_argument("--print-primes", action="store_true", help="imprimer la liste complète des m premiers (éviter pour y=79)")
    args = ap.parse_args()

    y = int(args.y)
    keep_primes = bool(args.print_primes or args.save_json)

    # borne simple utile (non imprimée si tu n'en as pas besoin)
    m_lb = int((y - 1) * math.log(2))
    print(f"\n=== calculate_m_requis ===")
    print(f"y = {y}  (premiers impairs ≥ y)")
    print(f"Borne rapide : m_requis(y) > (y-1)*ln(2) ≈ {(y-1)*math.log(2):.6f}  -> m > {m_lb}")

    res = compute_m_requis_and_certificate(y, keep_primes=keep_primes)
    print_compact_summary(res, reduced=args.reduced)

    if args.save_json:
        save_full_json(args.save_json, res, reduced=args.reduced)

if __name__ == "__main__":
    main()
