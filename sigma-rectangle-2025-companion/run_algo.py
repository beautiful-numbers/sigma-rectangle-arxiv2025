# run_algo_light.py
import math
from dataclasses import dataclass

# Échantillons pour calibrer kappa = m_min / (X^2/(2 ln X)) à partir de tes sorties
CALIB = [
    # (X, m_min_S_mesuré)
    (31,   238),
    (59,   591),
    (127,  2163),
    (151,  2977),
    (401,  16197),
    (601,  33025),
]

def kappa_from_samples(samples=CALIB):
    ratios = []
    for X, mmin in samples:
        val = (X*X)/(2.0*math.log(X))
        ratios.append(mmin / val)
    # médiane robuste
    ratios.sort()
    mid = len(ratios)//2
    return ratios[mid] if len(ratios)%2 else 0.5*(ratios[mid-1]+ratios[mid])

# kappa ~ 1.2–1.35 selon X ; on prend une médiane empirique.
KAPPA = kappa_from_samples()

def m_budget(D, X):
    # budget squarefree, nombre max de facteurs ≥ X à D chiffres
    return int((D*math.log(10.0)) / math.log(X))

def m_min_S_approx(X, kappa=KAPPA):
    # besoin minimal de facteurs pour franchir S>2 si p_min ≥ X (approx Mertens)
    return int(math.ceil(kappa * (X*X) / (2.0*math.log(X))))

@dataclass
class Crossing:
    D: int
    X_star: int
    kappa: float

def find_X_star(D, X_lo=11, X_hi=5000, kappa=KAPPA):
    # plus grand X tel que m_budget(D,X) >= m_min_S_approx(X)
    X_star = None
    for X in range(X_lo, X_hi+1):
        if m_budget(D,X) >= m_min_S_approx(X, kappa):
            X_star = X
    return X_star

def window_table(D, around=None, width=8, kappa=KAPPA):
    if around is None:
        around = find_X_star(D, kappa=kappa)
    if around is None:
        around = 200  # fallback
    rows = []
    for X in range(max(11, around-width), around+width+1):
        rows.append({
            "X": X,
            "m_budget": m_budget(D, X),
            "m_min_S_approx": m_min_S_approx(X, kappa),
            "status": "abondance possible" if m_budget(D,X) >= m_min_S_approx(X,kappa) else "déficience garantie"
        })
    return around, rows

if __name__ == "__main__":
    for D in [10_000, 50_000, 100_000, 500_000]:
        Xs = find_X_star(D)
        print(f"\n=== D = {D} ===")
        print(f"kappa (calibré) ≈ {KAPPA:.3f}")
        if Xs is None:
            print("Pas de croisement trouvé dans la fenêtre (augmente X_hi).")
            continue
        print(f"X*(D) ≈ {Xs}")
        Xs, rows = window_table(D, around=Xs)
        for r in rows:
            print(r)
