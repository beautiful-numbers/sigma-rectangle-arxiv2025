#!/usr/bin/env python3
# expand_until_pivot.py
#
# Compare "budget" vs "cost": on simule l'expansion contrainte d'un NPI
# via un graphe de dépendances r | sigma(p^e) et on mesure D_structure.
#
# Modes:
#   - factor  : calcule sigma(p^e), essaie d'en extraire des facteurs (trial division bornée),
#               et prend le cofacteur restant (souvent énorme) comme un "prime-like" (MR).
#   - zsig-lb : génère des "primes" forcés par Zsigmondy (candidats r ≡ 1 mod (e+1)), puis
#               promotion au prochain premier (next_prime). C'est une LB structurelle rapide.
#
# Stratégies d’exposants:
#   - square   : e = 2 pour tout p (optionnellement e=2 pour q si --force-e2-for-q)
#   - adaptive : e grossit jusqu’à e_max pour essayer de "sortir" plus de nouveaux facteurs
#
# Important:
#   - respect_y_rude : on n’ajoute que des r >= y_start.
#   - sample_frontier: traite à chaque depth seulement les X plus petits de la frontière pour tenir la charge.
#   - cap_sigma_digits: si digits(sigma(p^e)) > cap, on saute le factoring (pour rester fluide).
#
# Sortie: tableau "depth, omega, D_struct_min, D_struct_obs, newest_prime, status".
#
# Exemples:
#   python expand_until_pivot.py --y 79 --scenario B --exp-strategy square \
#       --respect-y-rude --m-requis 912 --max-depth 32 --print-header
#
#   python expand_until_pivot.py --y 101 --scenario A --mode zsig-lb \
#       --m-requis 1331 --max-depth 64 --zsig-per-node 2 --respect-y-rude --print-header

import argparse, math, sys, random

# ---------- util print ----------

def group_int(n: int) -> str:
    s = f"{n:,}".replace(",", " ")
    return s

# ---------- primalité & suivants ----------

_MR_BASES_64 = (2, 3, 5, 7, 11, 13, 17)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2,3,5,7,11,13,17,19,23,29)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    # Miller–Rabin
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    # pour > 2^64 on reste probabiliste, ça suffit pour notre usage
    bases = _MR_BASES_64 if n < (1<<64) else (2,3,5,7,11,13,17,19,23,29,31,37)
    for a in bases:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        ok = False
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                ok = True
                break
        if not ok:
            return False
    return True

def next_prime(n: int) -> int:
    if n <= 2:
        return 2
    p = n + 1 if (n % 2 == 0) else n
    while True:
        if is_probable_prime(p):
            return p
        p += 2

def next_prime_ge(n: int) -> int:
    if n <= 2:
        return 2
    p = n if n % 2 == 1 else n + 1
    while True:
        if is_probable_prime(p):
            return p
        p += 2

def first_prime_ge_1mod4(y: int) -> int:
    p = next_prime_ge(y)
    while p % 4 != 1:
        p = next_prime(p + 1)
    return p

# ---------- petits premiers pour factorisation bornée ----------

def simple_sieve(limit: int):
    bs = bytearray(b"\x01") * (limit + 1)
    bs[:2] = b"\x00\x00"
    r = int(limit**0.5)
    for p in range(2, r + 1):
        if bs[p]:
            step = p
            start = p * p
            bs[start:limit+1:step] = b"\x00" * (((limit - start)//step) + 1)
    return [i for i in range(2, limit + 1) if bs[i]]

_SMALL_PRIMES = simple_sieve(200000)  # suffisant pour "gratter" de petits facteurs

# ---------- sigma & outils ----------

def sigma_p_power(p: int, e: int) -> int:
    # sigma(p^e) = (p^(e+1)-1)//(p-1)
    return (pow(p, e + 1) - 1) // (p - 1)

def digits_est_sigma(p: int, e: int) -> int:
    # digits( (p^(e+1)-1)/(p-1) ) ~ floor((e+1) log10 p - log10(p-1)) + 1
    val = (e + 1) * math.log10(p) - math.log10(p - 1)
    return max(1, int(val) + 1)

def factor_smallish(n: int):
    """Trial division bornée par _SMALL_PRIMES. Renvoie une liste (multiset non nécessaire ici)."""
    out = []
    x = n
    for p in _SMALL_PRIMES:
        if p * p > x:
            break
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
    if x > 1:
        out.append(x)
    return out

# ---------- D_structure (observé & min) ----------

def D_struct_min_lb(y: int, omega: int) -> int:
    # ceil(omega * log10(y))
    if omega <= 0:
        return 0
    return int(math.ceil(omega * math.log10(y)))

def D_struct_obs_logsum(primes_set) -> int:
    """
    Robuste: ne compte que des entiers >= 2 (évite bool/None/0/1).
    Retourne 0 si vide; sinon floor(sum log10) + 1.
    """
    clean = []
    for p in primes_set:
        if isinstance(p, bool) or p is None:
            continue
        try:
            q = int(p)
        except Exception:
            continue
        if q >= 2:
            clean.append(q)
    if not clean:
        return 0
    s = 0.0
    for q in clean:
        # q>=2 garanti
        s += math.log10(q)
    return max(1, int(s) + 1)

# ---------- helpers candidats ----------

def ensure_valid_prime_like(r, y_min, promote_to_prime=True):
    """
    Normalise un 'candidat' r:
      - cast int,
      - impose r >= 2,
      - impose r >= y_min si respect-y-rude (le call le gère),
      - promotion au prochain premier si demandé.
    Retourne int valide ou None.
    """
    if r is None:
        return None
    try:
        r = int(r)
    except Exception:
        return None
    if r < 2:
        return None
    if promote_to_prime:
        r = next_prime_ge(r)
    return r

def zsig_lb_new_primes(p: int, e: int, count: int, y_min: int, respect_y_rude: bool) -> list:
    """
    Génère 'count' nouveaux 'primes' forcés par l’idée Zsigmondy:
    r ≡ 1 (mod e+1), r | p^(e+1)-1. On prend quelques constructions
    simples puis promotion au prochain premier. C’est une LB/heuristique.
    """
    k = e + 1
    cands = []
    # quelques formes croissantes & variées pour diversifier
    bases = [
        k * p + 1,
        k * p * p + 1,
        k * (p - 1) + 1,
        k * (p + 1) + 1,
        (k << 1) * p + 1,
        k * p * (k + 1) + 1,
    ]
    # boucle pour générer au moins 'count' candidats
    i = 0
    while len(cands) < count:
        r0 = bases[i % len(bases)]
        i += 1
        if respect_y_rude and r0 < y_min:
            r0 = y_min
        r = ensure_valid_prime_like(r0, y_min, promote_to_prime=True)
        if r is not None:
            cands.append(r)
    return cands[:count]

# ---------- expansion step ----------

def expand_step_factor(F: set, frontier: set, y_start: int, respect_y_rude: bool,
                       exp_strategy: str, force_e2_for_q: bool, e_max: int,
                       cap_sigma_digits: int | None, sample_frontier: int | None):
    """
    Mode factor: on calcule sigma(p^e), on tente de factoriser un peu,
    on pousse les cofacteurs comme nouveaux 'primes'.
    """
    if not frontier:
        return set(), None

    chosen = sorted(frontier)[: (sample_frontier or len(frontier))]
    next_frontier = set()
    newest = None

    for p in chosen:
        # Choix exposant
        if exp_strategy == "square":
            e = 2
            if force_e2_for_q:
                e = 2
        else:  # adaptive
            e = min(2, e_max)
        # on peut itérer un peu si adaptive
        e_list = [e] if exp_strategy == "square" else list(range(2, e_max + 1))

        for ee in e_list:
            # Cap sur digits(sigma) pour éviter des énormes entiers
            if cap_sigma_digits is not None:
                if digits_est_sigma(p, ee) > cap_sigma_digits:
                    continue
            try:
                sig = sigma_p_power(p, ee)
            except Exception:
                continue
            # facteurs (petits) + cofacteur
            facs = factor_smallish(sig)
            for r0 in facs:
                r = r0
                if respect_y_rude and r < y_start:
                    r = y_start
                r = ensure_valid_prime_like(r, y_start, promote_to_prime=True)
                if r is None:
                    continue
                if r not in F:
                    F.add(r)
                    next_frontier.add(r)
                    newest = r if (newest is None or r > newest) else newest
    return next_frontier, newest

def expand_step_zsig_lb(F: set, frontier: set, y_start: int, respect_y_rude: bool,
                        exp_strategy: str, e_max: int, zsig_per_node: int):
    """
    Mode zsig-lb: pas de factorisation, on force des nouveaux via zsig_lb_new_primes.
    """
    if not frontier:
        return set(), None

    chosen = sorted(frontier)  # déterministe
    next_frontier = set()
    newest = None

    # exposant utilisé: 2 si square; sinon on peut varier un peu
    e_list = [2] if exp_strategy == "square" else list(range(2, e_max + 1))

    for p in chosen:
        for ee in e_list:
            news = zsig_lb_new_primes(p, ee, zsig_per_node, y_start, respect_y_rude)
            for r in news:
                if respect_y_rude and r < y_start:
                    r = next_prime_ge(y_start)
                if r not in F:
                    F.add(r)
                    next_frontier.add(r)
                    newest = r if (newest is None or r > newest) else newest
    return next_frontier, newest

# ---------- scénarios init ----------

def init_scenario_A(y_start: int):
    # minimal: {y_start}
    return {next_prime_ge(y_start)}

def init_scenario_B(y_start: int):
    # forme d’Euler simplifiée: q ≡ 1 mod 4 (≥ y), plus 2 plus petits ≥ y (différents de q)
    q = first_prime_ge_1mod4(y_start)
    p1 = next_prime_ge(y_start)
    if p1 == q:
        p1 = next_prime(q + 1)
    p2 = next_prime(p1 + 1)
    return {q, p1, p2}, q

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--y", type=int, required=True, help="y_start (p_min ≥ y)")
    ap.add_argument("--scenario", choices=["A", "B"], default="A")
    ap.add_argument("--mode", choices=["factor", "zsig-lb"], default="factor")
    ap.add_argument("--exp-strategy", choices=["square", "adaptive"], default="square")
    ap.add_argument("--respect-y-rude", action="store_true")
    ap.add_argument("--force-e2-for-q", action="store_true")
    ap.add_argument("--e-max", type=int, default=6)
    ap.add_argument("--max-depth", type=int, default=16)
    ap.add_argument("--max-omega", type=int, default=5000)
    ap.add_argument("--m-requis", type=int, help="si fourni, calcule D_pivot=ceil(m*log10(y))")
    ap.add_argument("--print-header", action="store_true")
    ap.add_argument("--cap-sigma-digits", type=int, help="mode factor: saute factor si digits(sigma)>cap")
    ap.add_argument("--sample-frontier", type=int, help="traite X plus petits de la frontière par depth")
    ap.add_argument("--zsig-per-node", type=int, default=2, help="mode zsig-lb: nouveaux par noeud et exposant")
    args = ap.parse_args()

    y = args.y
    # coût
    if args.m_requis:
        D_pivot = int(math.ceil(args.m_requis * math.log10(y)))
        pivot_label = f"D_pivot (m_requis={args.m_requis}) = {group_int(D_pivot)}"
    else:
        # Borne constante rapide (m_lb ~ floor((y-1)*ln2)) ⇒ D_lb
        m_lb = int(math.floor((y - 1) * math.log(2.0)))
        D_pivot = int(math.ceil(m_lb * math.log10(y)))
        pivot_label = f"D_pivot_lower (const) = {group_int(D_pivot)}"

    # init
    if args.scenario == "A":
        F = init_scenario_A(y)
        q = None
    else:
        F, q = init_scenario_B(y)

    frontier = set(F)

    if args.print_header:
        print("=== expand_until_pivot ===")
        print(f"y = {y}  | scenario = {args.scenario}  | mode = {args.mode}  | strategy = {args.exp_strategy}  | max_depth = {args.max_depth}")
        print(f"{pivot_label}")
        print(f"respect_y_rude={args.respect_y_rude}  force_e2_for_q={args.force_e2_for_q}  e_max={args.e_max}")
        if args.mode == "factor":
            print(f"cap_sigma_digits={args.cap_sigma_digits or 'None'}  sample_frontier={args.sample_frontier or 'None'}")
        if args.mode == "zsig-lb":
            print(f"zsig_per_node={args.zsig_per_node}")
        print()
        print("depth |  omega  |  D_struct_min  |  D_struct_obs  |  newest_prime  |  status")

    # boucle
    newest = max(F) if F else None
    for depth in range(0, args.max_depth + 1):
        omega = len(F)
        Dmin = D_struct_min_lb(y, omega)
        Dobs = D_struct_obs_logsum(F)

        status = "budget < cost" if Dobs < D_pivot else "budget ≥ cost (hit pivot)"
        newest_str = f"{newest}" if newest is not None else "-"

        print(f"{depth:5d} | {omega:6d} | {Dmin:<14d} | {Dobs:<13d} | {newest_str:<14} |  {status}")

        if Dobs >= D_pivot:
            break
        if omega >= args.max_omega:
            print("[STOP] max_omega atteint.")
            break

        # étape d’expansion (sauf après la dernière ligne affichée)
        if depth == args.max_depth:
            break

        if args.mode == "factor":
            next_frontier, newest = expand_step_factor(
                F, frontier, y, args.respect_y_rude,
                args.exp_strategy, args.force_e2_for_q, args.e_max,
                args.cap_sigma_digits, args.sample_frontier
            )
        else:
            next_frontier, newest = expand_step_zsig_lb(
                F, frontier, y, args.respect_y_rude,
                args.exp_strategy, args.e_max, args.zsig_per_node
            )

        frontier = next_frontier
        if not frontier:
            print("\n[STOP] Frontière vide (aucun nouveau premier admissible).")
            break

if __name__ == "__main__":
    # lever toute limite de conversion int->str si Python 3.11+
    if hasattr(sys, "set_int_max_str_digits"):
        try:
            sys.set_int_max_str_digits(0)
        except Exception:
            pass
    main()
