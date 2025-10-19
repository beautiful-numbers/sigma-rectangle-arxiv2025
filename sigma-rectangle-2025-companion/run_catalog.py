# run_catalog.py
# ============================================================
# 1) Génère UNIQUEMENT les 27 parents via la minimalité "relèvement unitaire"
# 2) Pour chaque parent, construit une ancre MET déterministe :
#       - exposants pairs minimaux e_p = 2 jusqu'à σ(N)/N > 2 (glouton)
#       - pas besoin de "casser le carré" : MET fonctionne aussi pour N carré
# 3) Calcule σ(N), Delta = σ(N) - 2N (entiers), verdict
# 4) Produit :
#       - affichage console
#       - JSON canonique + sha256 (appendix_C_supplement.json)
#       - certificats lisibles (certificats_MET_ancres.txt)
# ============================================================

import json
import hashlib
from fractions import Fraction
from typing import Dict, Any, Tuple, List

from proof_algorithms import (
    generate_minimal_supports,
    ub_of_support,
)

# ---------- arithmétique ----------
def sigma_local(p: int, e: int) -> int:
    return (pow(p, e + 1) - 1) // (p - 1)

def local_ratio(p: int, e: int) -> Fraction:
    return Fraction(pow(p, e + 1) - 1, (p - 1) * pow(p, e))

def sigma_from_exponents(exps: Dict[int, int]) -> int:
    s = 1
    for p, e in exps.items():
        s *= sigma_local(p, e)
    return s

def sigma_ratio_from_exps(exps: Dict[int, int]) -> Fraction:
    r = Fraction(1, 1)
    for p, e in exps.items():
        r *= local_ratio(p, e)
    return r

# ---------- ancre MET : exposants pairs minimaux ----------
def minimal_even_exponents_for_abundance(
    S: Tuple[int, ...],
    target: Fraction = Fraction(2, 1),
    max_even: int = 12
) -> Dict[int, int]:
    exps = {p: 2 for p in S}
    ratio = sigma_ratio_from_exps(exps)
    if ratio > target:
        return exps

    def bump_gain(p: int, e: int) -> Fraction:
        num = local_ratio(p, e + 2)
        den = local_ratio(p, e)
        return Fraction(num, den)

    while ratio <= target:
        best_p, best_g = None, Fraction(1, 1)
        for p, e in exps.items():
            if e + 2 > max_even:
                continue
            g = bump_gain(p, e)
            if g > best_g:
                best_g, best_p = g, p
        if best_p is None:
            raise RuntimeError(
                "Impossible de dépasser 2 avec max_even=%d; augmente max_even."
                % max_even
            )
        exps[best_p] += 2
        ratio *= best_g
    return exps

# ---------- JSON canonique + hash ----------
def canonical_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))

def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

# ---------- format certificat ----------
def format_certificate_MET(cert_id: str, support: Tuple[int, ...], exps: Dict[int, int],
                           N: int, sigmaN: int, delta: int, s_ratio: Fraction) -> str:
    verdict = "abondant_strict" if delta > 0 else ("parfait" if delta == 0 else "deficient_strict")
    lines = []
    lines.append("="*100)
    lines.append(f"Certificat MET — ID: {cert_id}")
    lines.append("-"*100)
    lines.append(f"Support S: {support}")
    lines.append("Ancre (exposants pairs minimaux; MET fonctionne aussi si N est carré) :")
    lines.append(f"  exponents: {exps}")
    lines.append(f"  sigma(N)/N (exact) = {s_ratio.numerator}/{s_ratio.denominator}  (~ {float(s_ratio):.12f})")
    lines.append(f"  N = {N}")
    lines.append("")
    lines.append("MET (calcul exact) :")
    lines.append(f"  sigma(N) = {sigmaN}")
    lines.append(f"  2N       = {2*N}")
    lines.append(f"  Delta    = sigma(N) - 2N = {delta}")
    lines.append("")
    lines.append(f"Verdict via MET: {verdict}")
    return "\n".join(lines)

# ---------- main ----------
def main():
    print("Recherche des parents (supports minimaux, taille 8, sans 3, {5,7} forcés)…")
    parents = generate_minimal_supports()
    n_par = len(parents)
    print(f"Trouvés : {n_par} parents.")
    if n_par != 27:
        print("[NOTE] Le nombre attendu est 27 (Appendice C).")
        print("       Le code n'impose pas 27 ; si différent, vérifie P_MAX/environnement.\n")

    # Affichage
    print("="*150)
    print("Parents — Certificats MET (ancres déterministes)")
    print("="*150)
    header = f"{'ID':<6} | {'Support':<52} | {'UB(S)':<16} | {'Delta (ancre)':>24} | Verdict"
    print(header)
    print("-"*len(header))

    rows: List[dict] = []
    txt_blocks: List[str] = []

    for i, S in enumerate(parents, start=1):
        cert_id = f"C.{i:02d}"
        exps = minimal_even_exponents_for_abundance(S)
        N = 1
        for p, e in exps.items():
            N *= pow(p, e)
        sigmaN = sigma_from_exponents(exps)
        delta = sigmaN - 2 * N
        verdict = "abondant_strict" if delta > 0 else ("parfait" if delta == 0 else "deficient_strict")
        s_ratio = sigma_ratio_from_exps(exps)
        ub_val = ub_of_support(S)

        print(f"{cert_id:<6} | {str(S):<52} | {float(ub_val):<16.12f} | {delta:>24} | {verdict}")

        rec = {
            "id": cert_id,
            "support": list(S),
            "certificate": "MET(anchor)",
            "exponents": {str(p): int(e) for p, e in sorted(exps.items())},
            "N": str(N),
            "sigma": str(sigmaN),
            "Delta": str(delta),
            "S_ratio": {
                "numerator": str(s_ratio.numerator),
                "denominator": str(s_ratio.denominator)
            },
            "UB_support": {
                "numerator": str(ub_val.numerator),
                "denominator": str(ub_val.denominator)
            }
        }
        canon = canonical_json(rec)
        rec["sha256"] = sha256_hex(canon)
        rows.append(rec)

        txt_blocks.append(
            format_certificate_MET(cert_id, S, exps, N, sigmaN, delta, s_ratio)
        )

    # JSON
    json_fname = "appendix_C_supplement.json"
    with open(json_fname, "w", encoding="utf-8") as f:
        json.dump({"AppendixC_Parents_MET": rows}, f, sort_keys=True, ensure_ascii=False, indent=2)
    print(f"\nJSON écrit : {json_fname}")

    # Certificats texte
    txt_fname = "certificats_MET_ancres.txt"
    with open(txt_fname, "w", encoding="utf-8") as f:
        for block in txt_blocks:
            f.write(block + "\n")
    print(f"Certificats MET écrits : {txt_fname}")

if __name__ == "__main__":
    main()
