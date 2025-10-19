# run_ceiling_fast.py
# Objectif : déterminer rapidement si, à taille x et seuil premier y donnés,
# le "ceiling" d'abondance est franchissable.
# Deux modes O(1):
#   - const : test (y/(y-1))^m < 2 ?
#   - refined-approx : estimation via li-inverse et somme ~ log log(u) - log log(y)
#
# Remarque : ici on cherche un outil empirique rapide. Pas de produits entiers géants,
# pas d'énumération de millions de premiers. Tout est en logs, donc stable et rapide.

import math
import argparse
import re

LOG10_2 = math.log10(2.0)
LN_2 = math.log(2.0)
LN_10 = math.log(10.0)

def parse_x_or_digits(x_str: str | None, digits: int | None) -> int:
    """Retourne D = digits(x). Accepte --x sous forme '10^D' ou un entier,
    ou directement --digits D. Évite toute construction d'entier géant."""
    if digits is not None:
        if digits < 1:
            raise ValueError("digits doit être >= 1")
        return digits
    if x_str is None:
        raise ValueError("Spécifie --x ou --digits")
    m = re.fullmatch(r"\s*10\^(\d+)\s*", x_str)
    if m:
        return int(m.group(1)) + 1  # digits(10^D) = D+1
    # sinon on tente de compter les digits d'un entier raisonnable
    # (mais déconseillé pour les x énormes)
    try:
        x_val = int(x_str)
        if x_val <= 0:
            raise ValueError
        return len(str(x_val))
    except Exception:
        raise ValueError("Format --x non supporté (utilise '10^D' ou --digits D).")

def m_feasible(D: int, y: int) -> int:
    """m = floor(log_y x) avec x=10^(D-1..D). On prend x ≈ 10^(D-1) pour être sûr côté budget.
    Ici on adopte x = 10^(D-1) => ln x = (D-1) ln 10. C'est conservateur (légèrement plus petit)."""
    if y <= 2:
        raise ValueError("y doit être >= 3")
    ln_x = (D - 1) * LN_10
    return max(0, int(math.floor(ln_x / math.log(y))))

def const_decision(m: int, y: int):
    """Décide si (y/(y-1))^m < 2. On travaille en logs, donc O(1)."""
    if m == 0:
        return True, -LOG10_2  # produit = 1 < 2
    # log10((y/(y-1))^m / 2) = m*log10(1 + 1/(y-1)) - log10(2)
    # utiliser log1p pour la stabilité
    log10_ratio = m * (math.log1p(1.0 / (y - 1)) / math.log(10.0)) - LOG10_2
    return (log10_ratio < 0.0), log10_ratio

# ==== raffinage sans énumérer de premiers ====
# Idée : si on empile m premiers >= y, la somme des contributions
#   sum_{j=1..m} log(1 + 1/(p_j - 1))  ~  log log(u) - log log(y)
# où u ~ borne supérieure du m-ième premier au-delà de y, obtenu en inversant li.
# C'est heuristique mais très informatif pour de gros paramètres.
# On compare cette somme à ln 2.

def li_approx(u: float) -> float:
    """Approximation rapide de li(u) ~ u/ln u * (1 + 1/ln u + 2/ln^2 u).
    Suffisante pour le guidage ici (u grand)."""
    if u <= 3.0:
        # fallback grossier pour tout petit u
        return u / math.log(max(u, 2.001))
    L = math.log(u)
    invL = 1.0 / L
    return (u * invL) * (1.0 + invL + 2.0 * invL * invL)

def invert_li_increment(y: float, m: float) -> float:
    """Résout li(u) - li(y) = m pour u, par Newton amorti.
    On part d'un guess u0 ≈ y + m * ln y, puis on affine.
    Derivée : d/du li(u) = 1 / ln u."""
    if m <= 0:
        return y
    L_y = math.log(y)
    u = y + max(m * L_y, 10.0)  # guess
    for _ in range(8):
        Lu = math.log(u)
        f = li_approx(u) - li_approx(y) - m
        df = 1.0 / Lu
        step = f / df
        # amortissement pour robustesse
        u_next = u - 0.5 * step
        if u_next <= y + 1:
            u_next = 1.1 * (y + 1)
        # convergence grossière suffit
        if abs(u_next - u) / u < 1e-6:
            u = u_next
            break
        u = u_next
    return max(u, y + 1.0)

def refined_approx_decision(D: int, y: int, eps: float = 0.01):
    """Verdict heuristique sans premiers :
       - calcule m
       - approxime u tel que il y ait ~m premiers entre y et u
       - estime S ≈ log log(u) − log log(y)
       - compare S à ln 2 avec marge eps
    Retourne (verdict, delta, m, u)
      verdict ∈ {'likely_NO', 'likely_YES', 'uncertain'}
      delta = S - ln 2 (log naturel)
    """
    m = m_feasible(D, y)
    if m == 0:
        return 'likely_NO', -LN_2, 0, float(y)

    # approx u avec inversion rapide de li
    u = invert_li_increment(float(y), float(m))

    # estimation de la somme des contributions (log naturel)
    # prendre (u-1) et (y-1) pour coller à 1/(p-1)
    S_est = math.log(math.log(max(u - 1.0, 2.0))) - math.log(math.log(max(y - 1.0, 2.0)))
    delta = S_est - LN_2

    if delta < -eps:
        verdict = 'likely_NO'
    elif delta > eps:
        verdict = 'likely_YES'
    else:
        verdict = 'uncertain'
    return verdict, delta, m, u

def main():
    ap = argparse.ArgumentParser(description="Explorateur rapide du ceiling d'abondance (sans explosion de complexité).")
    ap.add_argument("--x", type=str, help="taille x (ex: 10^40000000).")
    ap.add_argument("--digits", type=int, help="ou bien directement digits(x).")
    ap.add_argument("--y", type=int, required=True, help="seuil des premiers (>=3).")
    ap.add_argument("--refined-approx", action="store_true", help="active l'estimation raffinée sans énumérer de premiers.")
    ap.add_argument("--eps", type=float, default=0.01, help="marge sur les logs naturels pour le verdict raffiné.")
    args = ap.parse_args()

    D = parse_x_or_digits(args.x, args.digits)
    y = args.y
    if y < 3:
        raise SystemExit("y doit être >= 3")

    print("\n=== Abundance ceiling (fast) ===")
    if args.x:
        print(f"x = {args.x}  |  digits(x) ≈ {D}")
    else:
        print(f"digits(x) = {D}")
    print(f"y = {y}")

    m = m_feasible(D, y)
    print("\n[Constant UB]")
    ok_const, log10_ratio = const_decision(m, y)
    print(f"  m = {m}")
    print(f"  décision: (y/(y-1))^m < 2 ?  {'YES' if ok_const else 'NO'}")
    print(f"  log10( (y/(y-1))^m / 2 ) ≈ {log10_ratio:.6f}  ({'< 0 => UB_const < 2' if log10_ratio < 0 else '≥ 0 => UB_const ≥ 2'})")

    if args.refined_approx:
        print("\n[Refined UB (approx, sans premiers)]")
        verdict, delta, m2, u = refined_approx_decision(D, y, eps=args.eps)
        # delta est en ln; afficher aussi en log10 pour intuition
        print(f"  m (budget) = {m2}")
        print(f"  borne u ~ {u:.3f}")
        print(f"  somme_est ≈ log(log(u-1)) - log(log(y-1))")
        print(f"  delta = somme_est - ln(2) ≈ {delta:.6f}  ({'<' if delta < 0 else '>'} 0)")
        print(f"  verdict ≈ {verdict}")

if __name__ == "__main__":
    main()
