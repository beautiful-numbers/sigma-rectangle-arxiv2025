# proof_algorithms.py
# ============================================================
# Objet : Générer UNIQUEMENT les 27 supports "parents" (taille 8, sans 3,
#         {5,7} forcés) tels que :
#   (i)   UB(T) >= 2
#   (ii)  pour tout p in T \ {5,7}, si q = next admissible prime > p
#         (impair, !=3, et q not in T), alors
#              UB( (T \ {p}) U {q} ) < 2
# => Minimalité par "relèvement unitaire" (anti-enfants).
#
# DFS avec élagage sûr :
#   - on complète pessimiste “au mieux” avec les plus petits primes restants ;
#     si même ce meilleur cas < 2 -> on coupe la branche.
# Caching agressif sur UB.
# ============================================================

from fractions import Fraction
from functools import lru_cache
from typing import List, Tuple

P_MAX = 20000  # large et sûr ; le DFS réel reste très petit car UB force des petits

# ---------- primalité simple (suffisant jusqu'à ~20000) ----------
def is_prime(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    if n % 3 == 0: return n == 3
    f, step = 5, 2
    while f * f <= n:
        if n % f == 0: return False
        f += step
        step = 6 - step
    return True

def next_prime(n: int) -> int:
    p = n + 1
    while not is_prime(p):
        p += 1
    return p

# ---------- construction des listes de premiers admissibles ----------
def admissible_primes() -> List[int]:
    # impairs != 3
    out = []
    for x in range(5, P_MAX + 1, 2):
        if x == 3:
            continue
        if is_prime(x):
            out.append(x)
    return out  # 5,7,11,13,17,...

# ---------- UB utilitaires ----------
def ub_ratio(p: int) -> Fraction:
    return Fraction(p, p - 1)

@lru_cache(maxsize=None)
def ub_of_support_cached(S_tuple: Tuple[int, ...]) -> Fraction:
    u = Fraction(1, 1)
    for p in S_tuple:
        u *= ub_ratio(p)
    return u

def ub_of_support(S: Tuple[int, ...]) -> Fraction:
    return ub_of_support_cached(tuple(sorted(S)))

def best_possible_completion_UB(current: Tuple[int, ...], tail: List[int], start_idx: int, need: int) -> Fraction:
    """
    Upper bound sûr : multiplie UB(current) par les 'need' plus petits
    premiers admissibles restants (>= tail[start_idx]) non déjà dans current.
    Si même ce "meilleur cas" est < 2 -> aucun achèvement n'atteindra 2.
    """
    u = ub_of_support(current)
    if need <= 0:
        return u
    cur_set = set(current)
    k = 0
    for i in range(start_idx, len(tail)):
        p = tail[i]
        if p in cur_set:
            continue
        u *= ub_ratio(p)
        k += 1
        if k == need:
            break
    # si on n'a pas réussi à prendre 'need' éléments, c'est mort de toute façon
    return u

# ---------- minimalité par relèvement unitaire ----------
def next_admissible_after(p: int, T_set: set) -> int:
    """
    Premier admissible q > p, impair, !=3, et q not in T_set.
    """
    q = p
    while True:
        q = next_prime(q)
        if q != 3 and (q % 2 == 1) and (q not in T_set):
            return q

def is_parent_minimal(T: Tuple[int, ...]) -> bool:
    """
    Critère parent-minimal (anti-enfants) :
      UB(T) >= 2
      et pour tout p in T\{5,7} :
         q = next admissible > p (q not in T)
         UB( (T\{p}) U {q} ) < 2
    """
    if ub_of_support(T) < 2:
        return False
    T_set = set(T)
    for p in T:
        if p in (5, 7):
            continue
        q = next_admissible_after(p, T_set)
        T2 = tuple(sorted((T_set - {p}) | {q}))
        if ub_of_support(T2) >= 2:
            return False
    return True

# ---------- génération (DFS) ----------
def generate_minimal_supports() -> List[Tuple[int, ...]]:
    ps = admissible_primes()
    if not (5 in ps and 7 in ps):
        raise RuntimeError("Premiers admissibles incomplets : 5 ou 7 manquent. Augmente P_MAX.")
    tail = [p for p in ps if p > 7]  # candidats après 7

    results: List[Tuple[int, ...]] = []
    seen = set()

    def dfs(current: Tuple[int, ...], start_idx: int):
        m = len(current)
        if m < 8:
            need = 8 - m
            # pruning "meilleur cas" (ajoute les 'need' plus petits restants)
            if best_possible_completion_UB(current, tail, start_idx, need) < 2:
                return
            for i in range(start_idx, len(tail)):
                p = tail[i]
                # croissant, pas de doublons
                if current and p <= current[-1]:
                    continue
                dfs(current + (p,), i + 1)
            return

        # m == 8
        if is_parent_minimal(current):
            key = tuple(sorted(current))
            if key not in seen:
                seen.add(key)
                results.append(key)

    dfs((5, 7), 0)
    results.sort()
    return results
