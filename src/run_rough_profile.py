# run_rough_profile.py
# ------------------------------------------------------------
# Portrait "y-rough" (aucun facteur premier ≤ y) pour x = 10^10 et 10^30,
# en retirant tous les premiers ≤ 999_983.
#
# Ce script :
#  1) estime Φ(x, y) = # { n ≤ x : p_min(n) > y } via la fonction de Buchstab ω(u),
#     avec u = log x / log y et l'approx classique Φ(x, y) ~ x * e^γ * ω(u) / log y ;
#  2) estime la part des premiers (π(x) − π(y)) via une approximation li(x) ;
#  3) estime la part des semi-premiers autorisés (p·q, y < p ≤ q, p·q ≤ x)
#     par intégration continue (densités ~ 1/log p, 1/log q) :
#        S ≈ ∫_{p=y}^{√x} [ Li(x/p) − Li(p) ] / log p  dp ;
#     (intégration en variable t = log p pour la stabilité numérique) ;
#  4) déduit la "profondeur combinatoire restante" :
#        m_max = ⌊ log x / log(y+1) ⌋
#     et la masse "autres" = Φ − 1 − premiers − semi-premiers (tronquée à 0).
#
# Conçu pour tourner sur une machine de bureau (ex: Intel i9 10th gen).
# Pas de dépendances externes.
# ------------------------------------------------------------

import math
from dataclasses import dataclass

# Constantes
EULER_GAMMA = 0.5772156649015328606065120900824024310421  # constante d'Euler–Mascheroni
EXP_EULER_GAMMA = math.e ** EULER_GAMMA                   # ≈ 1.78107

# ------------------------------------------------------------
# Approximation de li(x) et pi(x)
# ------------------------------------------------------------

def li_approx(x: float) -> float:
    """
    Approximation de la fonction logarithmique intégrale li(x) pour x >= 3,
    en utilisant l'expansion asymptotique :
      li(x) ~ x/log x * (1 + 1/log x + 2!/log^2 x + 3!/log^3 x + 4!/log^4 x)
    Suffisamment précise pour des ordres de grandeur (x jusqu'à 1e30).
    """
    if x < 3.0:
        # li(2) ~ 1.045... ; on donne une petite interpolation pour stabilité.
        # Mais dans nos usages, on ne s'en sert pratiquement pas pour x<3.
        return 1.04516378011749278484458888919 * (x - 1.0) / (3.0 - 1.0)
    L = math.log(x)
    invL = 1.0 / L
    # Somme des k! / L^k pour k=0..4
    s = 1.0 + invL + 2.0*invL**2 + 6.0*invL**3 + 24.0*invL**4
    return (x / L) * s

def pi_approx(x: float) -> float:
    """
    Approximation de π(x) par li(x) - li(2).
    Pour x très grand, la différence entre li et π est d'ordre subdominant.
    """
    if x < 2.0:
        return 0.0
    return max(li_approx(x) - li_approx(2.0), 0.0)

# ------------------------------------------------------------
# Buchstab ω(u) pour les y-rough : Φ(x, y) ~ x * e^γ * ω(u) / log y
# Définition :
#   pour 1 ≤ u ≤ 2, ω(u) = 1/u ;
#   pour u > 2, u ω'(u) + ω(u-1) = 0, avec continuité.
# On résout numériquement par pas d'Euler en u, suffisant ici (u ≤ ~5–8).
# ------------------------------------------------------------

def buchstab(u: float, du: float = 1e-3) -> float:
    """
    Approximation numérique de la fonction de Buchstab ω(u).
    Domaines utiles ici:
      - u = ln x / ln y pour (x, y) utilisés (ex: u≈1.666 pour 1e10/1e6, u≈5 pour 1e30/1e6).
    """
    if u < 1.0:
        return 0.0
    if u <= 2.0:
        return 1.0 / u

    # On résout ω sur [1, u] par pas 'du' :
    # Stockage sur grille : w[j] ≈ ω(uj), uj = 1 + j*du
    n = int(math.ceil((u - 1.0) / du))
    # Ajuster du pour que l'extrémité tombe pile
    du = (u - 1.0) / n if n > 0 else du

    # Initialisation sur [1,2] : ω(ξ) = 1/ξ
    # On stocke aussi pour ξ ∈ [1,2] car ω(u-1) y fera référence quand u ∈ [2,3].
    w = [0.0] * (n + 1)
    for j in range(n + 1):
        uj = 1.0 + j * du
        if uj <= 2.0:
            w[j] = 1.0 / uj
        else:
            break

    # Intégration pour u > 2 via u ω'(u) + ω(u-1) = 0 => ω'(u) = -ω(u-1)/u
    for j in range(1, n + 1):
        uj = 1.0 + j * du
        if uj <= 2.0:
            # déjà rempli par 1/uj
            continue
        # index correspondant à (u-1)
        um1 = uj - 1.0
        idx = int(round((um1 - 1.0) / du))
        idx = max(0, min(idx, len(w) - 1))
        w_um1 = w[idx]
        # pas d'Euler
        w[j] = w[j-1] + du * ( - w_um1 / uj )

    return w[-1]

def phi_rough_approx(x: float, y: float) -> float:
    """
    Approximation de Φ(x, y) = # { n ≤ x : p_min(n) > y }.
    Utilise la formule asymptotique avec Buchstab :
      Φ(x, y) ≈ x * e^γ * ω(u) / log y, u = log x / log y.
    """
    if x < 1.0:
        return 0.0
    if y < 2.0:
        # Rien n'est retiré : tous les n ≤ x (incluant 1)
        return x
    u = math.log(x) / math.log(y)
    w = buchstab(u)
    return (x * EXP_EULER_GAMMA * w) / math.log(y)

# ------------------------------------------------------------
# Estimation continue des semi-premiers autorisés (tous p, q > y, p ≤ q, pq ≤ x)
# S ≈ ∫_{p=y}^{√x} [ Li(x/p) − Li(p) ] / log p  dp
# On intègre en t = log p -> p = e^t, dp = e^t dt :
#   S ≈ ∫_{t=log y}^{log √x} [ Li(x/e^t) − Li(e^t) ] / t * e^t  dt
# Trapèzes sur une grille en t (log-équidistante en p).
# ------------------------------------------------------------

def semiprimes_truncated_estimate(x: float, y: float, grid: int = 4000) -> float:
    """
    Estime le nombre de semi-premiers ≤ x dont les deux facteurs premiers sont > y.
    Utilise une intégration continue avec densité PNT (1/log).
    """
    if y * y > x:
        return 0.0

    t_lo = math.log(y)
    t_hi = 0.5 * math.log(x)  # log(√x)

    # Grille en t
    grid = max(2, grid)
    dt = (t_hi - t_lo) / (grid - 1)
    total = 0.0

    for i in range(grid):
        t = t_lo + i * dt
        p = math.exp(t)
        # Integrand f(p) = [ Li(x/p) − Li(p) ] / log p
        # En variables t : f_t(t) = ([ Li(x/e^t) − Li(e^t) ] / t) * e^t
        if t <= 0.0:  # sécurité
            continue
        Li_upper = li_approx(x / p)
        Li_lower = li_approx(p)
        inner = max(Li_upper - Li_lower, 0.0)
        val = (inner / t) * p
        # Poids trapèze
        w = 0.5 if (i == 0 or i == grid - 1) else 1.0
        total += w * val

    return max(total * dt, 0.0)

# ------------------------------------------------------------
# Analyse complète pour un couple (x, y)
# ------------------------------------------------------------

@dataclass
class RoughPortrait:
    x: int
    y: int
    u: float
    buchstab_omega: float
    phi_rough: float
    primes_allowed: float
    semiprimes_allowed: float
    others_allowed: float
    m_max: int

def analyze_case(x: int, y: int, semiprime_grid: int = 4000) -> RoughPortrait:
    """
    Calcule les métriques “y-rough” ≤ x :
      - Φ(x, y) (approx)
      - #primes (π(x) − π(y)) (approx)
      - #semi-premiers autorisés (approx)
      - "autres" = Φ − 1 − primes − semi (tronqué à 0)
      - m_max = ⌊ log x / log(y+1) ⌋
    """
    x = int(x)
    y = int(y)

    # u et ω(u)
    u = math.log(x) / math.log(y)
    w = buchstab(u)

    # Φ(x, y)
    phi = phi_rough_approx(x, y)

    # Primes autorisés
    primes = max(pi_approx(x) - pi_approx(y), 0.0)

    # Semi-premiers autorisés
    semi = semiprimes_truncated_estimate(float(x), float(y), grid=semiprime_grid)

    # m_max : profondeur combinatoire possible (distincts > y)
    m_max = int(math.floor(math.log(x) / math.log(y + 1.0)))

    # "Autres" (≥ 3 facteurs) : tronqué à 0 (ceci reste un proxy)
    others = max(phi - 1.0 - primes - semi, 0.0)

    return RoughPortrait(
        x=x,
        y=y,
        u=u,
        buchstab_omega=w,
        phi_rough=phi,
        primes_allowed=primes,
        semiprimes_allowed=semi,
        others_allowed=others,
        m_max=m_max
    )

# ------------------------------------------------------------
# Affichage lisible
# ------------------------------------------------------------

def human(x: float) -> str:
    """ Formatage compact avec 3 sig. figs. """
    if x < 1e3:
        return f"{x:.3g}"
    mag = int(math.floor(math.log10(abs(x)))) if x != 0 else 0
    base = x / (10 ** mag)
    return f"{base:.3f}e{mag:+d}"

def print_portrait(p: RoughPortrait) -> None:
    print("\n------------------------------------------------------------")
    print(f"y-rough portrait  |  x = {p.x:.0f}  |  y = {p.y}")
    print("------------------------------------------------------------")
    print(f"u = log x / log y         : {p.u:.6f}")
    print(f"Buchstab ω(u)             : {p.buchstab_omega:.6f}")
    print(f"Φ(x, y) (y-rough, approx) : {human(p.phi_rough)}")
    print(f"  ├─ premiers autorisés   : {human(p.primes_allowed)}")
    print(f"  ├─ semi-premiers (>y)   : {human(p.semiprimes_allowed)}")
    print(f"  └─ autres (≥ 3 facteurs): {human(p.others_allowed)}")
    print(f"Profondeur m_max (distincts > y) : {p.m_max}")
    if p.y * p.y > p.x:
        print("Note : x < y^2 → aucun composé autorisé (seulement 1 et des premiers).")
    print("------------------------------------------------------------")

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    # Paramètre de retrait des petits premiers
    y = 999_983

    # Cas demandés
    xs = [10**10, 10**30]

    # Grille d'intégration pour les semi-premiers (compromis précision/temps)
    grid = 4000  # ~ quelques centaines de ms sur i9 10th gen

    print("Rough-combinatorics under prime removal")
    print(f"(all primes ≤ {y} removed; look at y-rough numbers ≤ x)")
    print("Approximations used: Buchstab ω(u) (Euler step), li-approx (PNT series),")
    print("                     continuous semiprime integral in log-space.")

    for x in xs:
        p = analyze_case(x, y, semiprime_grid=grid)
        print_portrait(p)

if __name__ == "__main__":
    main()
