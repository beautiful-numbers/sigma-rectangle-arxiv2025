#!/usr/bin/env python3
# run_seuil_abondance.py
# Explorer le "seuil d'abondance" fenêtré : pour un x donné, trouver y tel que
# tout n <= x avec p_min(n) >= y est déficient (UB affiné < 2).
#
# Décisions exactes en entiers (on compare ∏p et 2∏(p-1)).
# Aucune conversion d'entiers géants en chaîne : pas de limite "4300 digits".
#
# Exemples :
#   python run_seuil_abondance.py --x 10^3000 --y 79 --refined --print-primes
#   python run_seuil_abondance.py --x 10^4299 --find-refined-ceiling
#   python run_seuil_abondance.py --x 10^80 --y 1000000000019 --refined

import argparse, math, sys
from decimal import Decimal, getcontext

# lever toute limite de conversion si jamais on imprimait un entier
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

LOG10 = math.log(10.0)

# ---------- parsing de x et outils log ----------
def parse_x(expr: str):
    """
    renvoie un dict avec :
      - ln_x : ln(x) en double (suffisant pour m = floor(ln x / ln y))
      - digits_est : estimation des digits de x (pour affichage)
      - label : rendu court (ex "10^3000")
    accepte formats :
      "10^D"  |  "10**D"  |  entier décimal  |  float scientifique (usage léger)
    """
    s = expr.strip().lower().replace(" ", "")
    if s.startswith("10^"):
        D = int(s[3:])
        ln_x = D * LOG10
        digits_est = D + 1
        label = f"10^{D}"
        return dict(ln_x=ln_x, digits_est=digits_est, label=label)
    if s.startswith("10**"):
        D = int(s[4:])
        ln_x = D * LOG10
        digits_est = D + 1
        label = f"10^{D}"
        return dict(ln_x=ln_x, digits_est=digits_est, label=label)
    # entier décimal brut
    if s.isdigit():
        # pas de str(n) ici ; on calcule ln via Decimal si besoin
        getcontext().prec = max(50, len(s) // 3 + 20)
        ln10 = Decimal(10).ln()
        # ln(n) = ln( (premiers_k chiffres) * 10^(len-...)) — on passe direct par Decimal
        n_dec = Decimal(s)
        ln_x_dec = n_dec.ln()
        ln_x = float(ln_x_dec)  # conversion pour m-feasible
        # estimation digits
        digits_est = int(math.floor(ln_x / LOG10)) + 1
        label = expr
        return dict(ln_x=ln_x, digits_est=digits_est, label=label)
    # fallback: float
    x = float(expr)
    ln_x = math.log(x)
    digits_est = int(math.floor(ln_x / LOG10)) + 1
    return dict(ln_x=ln_x, digits_est=digits_est, label=expr)

def m_budget(ln_x: float, y: int) -> int:
    # m = floor( ln x / ln y )
    return int(math.floor(ln_x / math.log(y)))

# ---------- génération de premiers ----------
def simple_sieve(limit: int):
    bs = bytearray(b"\x01") * (limit + 1)
    bs[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if bs[p]:
            step = p
            start = p * p
            bs[start:limit+1:step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if bs[i]]

# Miller–Rabin déterministe pour n < 2^64
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
    # n-1 = d * 2^s
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    for a in _MR_BASES_64:
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
    génère les m premiers nombres premiers >= y.
    si y < 10^7, on sieve jusqu'à une borne raisonnable,
    sinon on avance avec Miller–Rabin.
    """
    out = []
    if y <= 10_000_000:
        # borne pratique : on surestime le m-ième premier après y
        # p_{k+m} ~ (y + m*(log y + log log y)) grossièrement
        ly = math.log(max(3, y))
        bound = int(y + m * (ly + math.log(ly + 1)))
        bound = max(bound, y + 2000)
        primes = simple_sieve(bound)
        # sauter ceux < y
        lo = 0
        # trouver premier index >= y
        import bisect
        lo = bisect.bisect_left(primes, y)
        out = primes[lo:lo + m]
        if len(out) < m:
            # fallback : compléter en MR
            last = primes[-1] if primes else y
            p = last + 1
            if p % 2 == 0:
                p += 1
            while len(out) < m:
                if is_probable_prime(p):
                    out.append(p)
                p += 2
        return out
    else:
        p = y if is_probable_prime(y) else next_prime(y)
        out.append(p)
        while len(out) < m:
            p = next_prime(p + 1)
            out.append(p)
        return out

# ---------- UB constant et UB raffiné ----------
def ub_const_lt2(ln_x: float, y: int):
    """
    décide avec la borne constante :
      m = floor(ln x / ln y)
      test (y/(y-1))^m < 2 ?
    renvoie (m, decision_bool, margin_log10)
      margin_log10 = log10( (y/(y-1))^m / 2 )
    """
    m = m_budget(ln_x, y)
    if m <= 0:
        return (m, True, -1.0)  # produit vide : 1 < 2
    ratio = y / (y - 1.0)
    # utiliser log10 pour afficher une marge lisible
    margin = m * math.log10(ratio) - math.log10(2.0)
    return (m, margin < 0.0, margin)

def ub_refined_cmp(y: int, m: int, want_primes=False):
    """
    calcule le produit exact ∏_{j=1..m} p_j/(p_j-1) avec p_j >= y
    renvoie (crossed, count, maybe_primes, digits_num_approx, digits_den_approx)
      crossed : True si produit >= 2, sinon False
      count   : nombre de facteurs réellement empilés (m si on va au bout)
      maybe_primes : la liste des p_j si want_primes
      digits_*_approx : ~floor(sum log10) + 1 (pas de str géant)
    comparaison exacte : on maintient num = ∏p, den = ∏(p-1) et on teste num >= 2*den
    """
    if m <= 0:
        return (False, 0, [] if want_primes else None, 1, 1)
    primes = primes_from_y(y, m)
    num = 1
    den = 1
    sumlog_num = 0.0
    sumlog_den = 0.0
    for idx, p in enumerate(primes, 1):
        num *= p
        den *= (p - 1)
        sumlog_num += math.log10(p)
        sumlog_den += math.log10(p - 1)
        # test exact d'arrêt
        if num >= 2 * den:
            dnum = int(sumlog_num) + 1
            dden = int(sumlog_den) + 1
            return (True, idx, primes if want_primes else None, dnum, dden)
    dnum = int(sumlog_num) + 1
    dden = int(sumlog_den) + 1
    return (False, m, primes if want_primes else None, dnum, dden)

# ---------- recherche du seuil raffiné ----------
def prime_iter(start: int):
    p = start if is_probable_prime(start) else next_prime(start)
    while True:
        yield p
        p = next_prime(p + 1)

def find_refined_ceiling(ln_x: float, y_hint: int = 5, y_max: int = 1_000_000):
    """
    renvoie le plus petit y tel que UB raffiné < 2
    stratégie :
      1) on part d’un y bas et on cherche une borne haute YH avec produit < 2
         (on augmente y par sauts exponentiels jusqu’à trouver un échec)
      2) binaire sur les nombres premiers entre basse et haute bornes
    """
    # helper: décide pour un y
    def decide(y: int):
        m = m_budget(ln_x, y)
        crossed, cnt, _, _, _ = ub_refined_cmp(y, m, want_primes=False)
        return (m, crossed)

    # étape 1 : trouver une borne haute où ça NE dépasse pas 2
    yl = max(3, y_hint)
    yh = yl
    m, crossed = decide(yh)
    steps = 0
    while crossed and yh < y_max and steps < 40:
        yh = next_prime(int(yh * 1.5) + 1)
        m, crossed = decide(yh)
        steps += 1
    if crossed:
        # pas trouvé de coupe, renvoyer None
        return None

    # maintenant il existe yl' < yh tel que à yl' ça dépasse, à yh non
    # on fixe yl comme dernier y qui dépasse (on recule par paliers)
    yl = 3
    # on monte yl jusqu’à trouver "crossed" puis on le garde
    it = prime_iter(yl)
    yl = next(it)  # premier >= 3
    last_cross = 3
    while yl < yh:
        m, crossed = decide(yl)
        if crossed:
            last_cross = yl
            yl = next(it)
        else:
            break
    low = last_cross
    high = yh

    # binaire sur les premiers entre low et high
    # on construit la liste des premiers dans l'intervalle pour un bisect manuel
    plist = []
    p = low
    while p <= high:
        plist.append(p)
        p = next_prime(p + 1)

    lo = 0
    hi = len(plist) - 1
    ans = None
    while lo <= hi:
        mid = (lo + hi) // 2
        y_try = plist[mid]
        m, crossed = decide(y_try)
        if crossed:
            lo = mid + 1
        else:
            ans = y_try
            hi = mid - 1
    return ans

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Abundance ceiling explorer (fenêtré)")
    ap.add_argument("--x", required=True, help="taille max (ex: 10^3000, 10**80, entier)")
    ap.add_argument("--y", type=int, help="seuil y (p_min >= y)")
    ap.add_argument("--refined", action="store_true", help="utiliser UB raffiné avec vrais premiers >= y")
    ap.add_argument("--print-primes", action="store_true", help="afficher la liste des premiers utilisés")
    ap.add_argument("--find-ceiling", action="store_true", help="trouver y* avec UB constant (rapide)")
    ap.add_argument("--find-refined-ceiling", action="store_true", help="trouver y_ref avec UB raffiné (binaire)")
    args = ap.parse_args()

    X = parse_x(args.x)
    ln_x = X["ln_x"]
    print("\n=== Abundance ceiling explorer ===")
    print(f"x = {X['label']}  |  digits(x) = {X['digits_est']}\n")

    if args.find_ceiling:
        # cherche y* minimal (version constante)
        # on balaye y croissant jusqu'à (y/(y-1))^m < 2
        y = 3
        y_star = None
        while y < 2_000_000:
            m, ok, margin = ub_const_lt2(ln_x, y)
            if ok:
                y_star = y
                break
            y = next_prime(y + 1)
        if y_star is None:
            print("[Ceiling-const] pas trouvé dans la fenêtre.")
        else:
            print("[Ceiling-const] plus petit y tel que UB_const < 2 :")
            print(f"  y_const = {y_star}")
            m, ok, margin = ub_const_lt2(ln_x, y_star)
            print(f"  m = floor(log_y x) = {m}")
            print(f"  log10( (y/(y-1))^m / 2 ) ≈ {margin:.6f}  (< 0 implique UB_const < 2)")
        return

    if args.find_refined_ceiling:
        y_ref = find_refined_ceiling(ln_x, y_hint=5, y_max=10_000_000)
        if y_ref is None:
            print("[Ceiling-refined] pas de y ≤ 10^7 forçant UB raffiné < 2 pour ce x.")
        else:
            m = m_budget(ln_x, y_ref)
            crossed, cnt, _, dnum, dden = ub_refined_cmp(y_ref, m, want_primes=False)
            # crossed doit être False (coupe)
            print("[Ceiling-refined] plus petit y tel que UB raffiné < 2 :")
            print(f"  y_ref = {y_ref}")
            print(f"  m = floor(log_y x) = {m}")
            print(f"  décision: produit < 2  (comparaison exacte ∏p vs 2∏(p-1))")
            print(f"  digits approx: ∏p ≈ {dnum},  2∏(p-1) ≈ {dden}")
        return

    if args.y is None:
        print("fournis --y ou utilise --find-ceiling / --find-refined-ceiling")
        return

    y = int(args.y)
    print("--- Parameters ---")
    print(f"y = {y}")
    m = m_budget(ln_x, y)
    print(f"m_feasible (exact) = floor(log_y x) = {m}\n")

    # UB constant (utile comme repère)
    m_const, ok_const, margin = ub_const_lt2(ln_x, y)
    print("[Constant UB]")
    print(f"  m = {m_const}")
    print(f"  décision: (y/(y-1))^m < 2 ?  {'YES' if ok_const else 'NO'}")
    print(f"  log10( (y/(y-1))^m / 2 ) ≈ {margin:.6f}  (< 0 implique UB_const < 2)\n")

    if args.refined:
        crossed, cnt, plist, dnum, dden = ub_refined_cmp(y, m, want_primes=args.print_primes)
        print("[Refined UB with real primes ≥ y]")
        print(f"  m (primes stacked) = {cnt}")
        if crossed:
            print("  décision: ∏ p/(p-1) ≥ 2 ?  YES")
            print("  check exact: ∏p ≥ 2∏(p-1)")
        else:
            print("  décision: ∏ p/(p-1) ≥ 2 ?  NO")
            print("  check exact: ∏p < 2∏(p-1)")
        print(f"  digits approx: ∏p ≈ {dnum},  2∏(p-1) ≈ {dden}")

if __name__ == "__main__":
    main()
