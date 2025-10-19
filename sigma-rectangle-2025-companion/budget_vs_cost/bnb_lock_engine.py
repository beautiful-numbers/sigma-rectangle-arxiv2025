#!/usr/bin/env python3
# bnb_lock_engine.py  — Version optimisée (Axes 1–3)
#
# Objet :
#   Annexes (4.8) – Verrous structurels et couverture effective.
#   Implémente les optimisations :
#     Axe 1  : Lazy & Cached (primes pré-calculés, factorisations/ordres/σ mémorisés,
#               UB_prefix en O(1), factorisation rapide via Pollard–Rho si nécessaire)
#     Axe 2  : Fail Fast (choix heuristique de la base, des exposants, et de r “balle empoisonnée”)
#     Axe 3  : Strongest Weapon First (ordre des coupes : v2 -> S×UB -> L1 -> L2 -> L3)
#
#   - M(y) : PPCM des ordres admis (k = e+1) jusqu’à K_max -> classes r ≡ 1 (mod d)
#   - D_ram(y) : # { r < y premier : ∃ d | M(y) avec r ≡ 1 (mod d) }
#   - C_par(y) : borne simple du # de bases avec e pair (profil Euler) = ω_cap - 1
#   - UB_m(y) (préfixe) et test S_part × UB_restant < 2 (UB affiné depuis X_local)
#   - Moteur BnB + coupes L :
#       L0 : v2(σ(n)) > 1  (coupe immédiate et très bon marché)
#       L1 : apparition d’un r < y (y-rude violé)
#       L2 : cycles interdits (inclut conflit Euler)
#       L3 : neutralité q-adique impossible (cadre générique)
#   - Exceptions de Zsigmondy traitées explicitement.
#
# Exemples :
#   1) Tableau M(y)/D_ram(y) pour y = 3,5,7,11 (K_max=32) :
#      python bnb_lock_engine.py table --ys 3,5,7,11 --kmax 32
#
#   2) UB et C_par :
#      python bnb_lock_engine.py table --ymin 3 --ymax 31 --kmax 40 --omega-cap 8
#
#   3) BnB+L (journal succinct) pour y = 11, ω_cap = 22 :
#      python bnb_lock_engine.py bnb --y 11 --omega-cap 22 --even-exps 2 --euler-exps 1 --max-nodes 5000 --sieve-max 200000
#
# NB :
#   - Pas d’accès réseau. Pollard–Rho + MR pour factoriser proprement quand nécessaire.
#   - Tous les flottants sont en double ; suffisant pour les UB.
#   - Pour des runs “démo rapides”, réduire les sets d’exposants (ex: --even-exps 2 --euler-exps 1).

import argparse
import math
import random
import sys
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set, Iterable
import csv
import bisect

# ============================================================
# ============   AXE 1 : Lazy & Cached (global)   ============
# ============================================================

# ---- Caches arithmétiques ----
_SIG_CACHE: Dict[Tuple[int,int], int] = {}
_FAC_CACHE: Dict[int, Dict[int,int]] = {}
_ORD_CACHE: Dict[Tuple[int,int], Optional[int]] = {}
_SFAC_CACHE: Dict[Tuple[int,int], float] = {}
_V2SIG_CACHE: Dict[Tuple[int,int], int] = {}
_VQ_CACHE: Dict[Tuple[int,int,int], int] = {}
_TAU_CACHE: Dict[int, int] = {}

# ---- Primalité (Miller–Rabin) & Pollard–Rho ----

def _mulmod(a: int, b: int, m: int) -> int:
    return (a * b) % m

def _powmod(a: int, d: int, m: int) -> int:
    return pow(a, d, m)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    # petits cas
    small_primes = (2,3,5,7,11,13,17,19,23,29)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
    # bases MR (suffisantes pour 64-bit ; pour grands n, proba négligeable)
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29):
        if a % n == 0:
            continue
        x = _powmod(a, d, n)
        if x == 1 or x == n - 1:
            continue
        skip = False
        for _ in range(s - 1):
            x = _mulmod(x, x, n)
            if x == n - 1:
                skip = True
                break
        if skip:
            continue
        return False
    return True

def _pollard_rho(n: int) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = random.randrange(1, n-1)
        f = lambda x: (x * x + c) % n
        x, y, d = 2, 2, 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d

def factorint(n: int, fac: Optional[Dict[int,int]] = None) -> Dict[int,int]:
    """Factorisation complète (MR + Pollard–Rho), avec cache."""
    if n in _FAC_CACHE:
        return dict(_FAC_CACHE[n])
    if fac is None:
        fac = {}
    if n <= 1:
        _FAC_CACHE[n] = dict(fac)
        return dict(fac)
    if is_probable_prime(n):
        fac[n] = fac.get(n, 0) + 1
        _FAC_CACHE[n] = dict(fac)
        return dict(fac)
    d = _pollard_rho(n)
    factorint(d, fac)
    factorint(n // d, fac)
    _FAC_CACHE[n] = dict(fac)
    return dict(fac)

# ---- Crible & gestion des premiers ----

def sieve_primes(limit: int) -> List[int]:
    if limit < 2:
        return []
    bs = bytearray(b"\x01")*(limit+1)
    bs[:2] = b"\x00\x00"
    r = int(limit**0.5)
    for p in range(2, r+1):
        if bs[p]:
            step = p
            start = p*p
            bs[start:limit+1:step] = b"\x00"*(((limit - start)//step)+1)
    return [i for i in range(2, limit+1) if bs[i]]

def extend_primes_from(last: int, count: int) -> List[int]:
    """Génère `count` prochains premiers impairs > last (fallback si X dépasse le crible)."""
    out = []
    n = last + 1 if last % 2 == 0 else last + 2
    while len(out) < count:
        if is_probable_prime(n):
            out.append(n)
        n += 2
    return out

# ============================================================
# =================  Outils arithmétiques  ===================
# ============================================================

def lcm(a: int, b: int) -> int:
    return abs(a*b)//math.gcd(a,b)

def lcm_many(vals: Iterable[int]) -> int:
    L = 1
    for v in vals:
        L = lcm(L, v)
    return L

def divisors(n: int) -> List[int]:
    fac = factorint(n)
    ds = [1]
    for p, e in fac.items():
        cur = []
        pe = 1
        for _ in range(e):
            pe *= p
            for d in ds:
                cur.append(d * pe)
        ds += cur
    ds = list(sorted(set(ds)))
    return ds

def tau_of(n: int) -> int:
    t = _TAU_CACHE.get(n)
    if t is not None:
        return t
    fac = factorint(n)
    v = 1
    for e in fac.values():
        v *= (e+1)
    _TAU_CACHE[n] = v
    return v

def sigma_p_pow_e(p: int, e: int) -> int:
    key = (p,e)
    s = _SIG_CACHE.get(key)
    if s is None:
        s = (pow(p, e+1) - 1)//(p-1)
        _SIG_CACHE[key] = s
    return s

def S_factor(p: int, e: int) -> float:
    key = (p,e)
    v = _SFAC_CACHE.get(key)
    if v is None:
        v = sigma_p_pow_e(p,e) / float(pow(p, e))
        _SFAC_CACHE[key] = v
    return v

def v2_of(n: int) -> int:
    return (n & -n).bit_length() - 1 if n else 0

def v2_sigma_p_e(p: int, e: int) -> int:
    """v2(σ(p^e)) via LTE (p impair). σ(p^e) = (p^{m}-1)/(p-1), m=e+1.
       - si m impair : v2 = v2(m)
       - si m pair   : v2 = v2(p+1) + v2(m) - 1
    """
    key = (p,e)
    got = _V2SIG_CACHE.get(key)
    if got is not None:
        return got
    m = e + 1
    if m % 2 == 1:
        v = v2_of(m)
    else:
        v = v2_of(p + 1) + v2_of(m) - 1
    if v < 0:
        v = 0
    _V2SIG_CACHE[key] = v
    return v

def order_mod(a: int, mod: int) -> Optional[int]:
    key = (a, mod)
    t = _ORD_CACHE.get(key)
    if t is not None:
        return t
    if math.gcd(a, mod) != 1:
        _ORD_CACHE[key] = None
        return None
    phi = mod - 1
    for d in divisors(phi):
        if pow(a, d, mod) == 1:
            _ORD_CACHE[key] = d
            return d
    _ORD_CACHE[key] = None
    return None

def v_q_of_sigma_p_e(q: int, p: int, e: int) -> int:
    key = (q,p,e)
    v = _VQ_CACHE.get(key)
    if v is not None:
        return v
    s = sigma_p_pow_e(p,e)
    c = 0
    while s % q == 0:
        s //= q
        c += 1
    _VQ_CACHE[key] = c
    return c

# ============================================================
# ===================  UB et gestion primes  =================
# ============================================================

@dataclass
class UBPrimes:
    start_y: int
    sieve_max: int
    primes_list: List[int] = field(default_factory=list)      # tous les premiers (impairs inclus)
    primes_ge_y: List[int] = field(default_factory=list)      # premiers impairs >= y
    ub_prefix: List[float] = field(default_factory=list)      # préfixe de ∏ p/(p-1) sur primes_ge_y

    def ensure_length_from(self, idx: int, need_more: int):
        """S'assure que primes_ge_y possède au moins idx+need_more éléments ; étend si besoin."""
        want = idx + need_more
        if want <= len(self.primes_ge_y):
            return
        # Étendre au-delà du crible avec MR si besoin
        last = self.primes_ge_y[-1] if self.primes_ge_y else max(3, self.start_y)
        extra = extend_primes_from(last, want - len(self.primes_ge_y))
        for p in extra:
            if p % 2 == 0:  # on garde les impairs
                continue
            self.primes_ge_y.append(p)
            # maj ub_prefix
            last_prod = self.ub_prefix[-1] if self.ub_prefix else 1.0
            self.ub_prefix.append(last_prod * (p/(p-1)))

    def ub_from_X(self, X: int, m: int) -> float:
        """UB_m mais à partir de X (≥ y) : ∏_{j=0..m-1} p_j/(p_j-1) pour les m plus petits p≥X."""
        if m <= 0:
            return 1.0
        # index binaire dans primes_ge_y
        if not self.primes_ge_y:
            self.ensure_length_from(0, 1)
        idx = bisect.bisect_left(self.primes_ge_y, max(self.start_y, X))
        self.ensure_length_from(idx, m)
        # produit préfixe : UB = prefix[idx+m] / prefix[idx]
        prefix = self.ub_prefix
        num = prefix[idx + m - 1]
        den = prefix[idx - 1] if idx > 0 else 1.0
        return num / den

# ============================================================
# ===================   M(y) / D_ram(y)   ====================
# ============================================================

def build_M_y(K_max: int, allow_even_exps=True, allow_one_odd_euler=True) -> int:
    Ks = set()
    if allow_even_exps:
        for e in range(2, K_max+1):
            if e % 2 == 0:
                Ks.add(e+1)
    if allow_one_odd_euler:
        for e in range(1, K_max+1, 2):
            if e % 4 == 1:
                Ks.add(e+1)
    Ks = [k for k in Ks if k >= 2 and k <= (K_max+1)]
    if not Ks:
        return 1
    return lcm_many(Ks)

def primes_up_to(n: int) -> List[int]:
    return sieve_primes(n)

def D_ram_of_y(y: int, M_y: int) -> Tuple[int, Dict[int,int]]:
    if y <= 3:
        return (0, {})
    r_primes = primes_up_to(y-1)
    Ds = [d for d in divisors(M_y) if d >= 2]
    cover_set: Set[int] = set()
    breakdown: Dict[int,int] = {}
    for d in Ds:
        c = 0
        for r in r_primes:
            if r % d == 1:
                cover_set.add(r)
                c += 1
        breakdown[d] = c
    return (len(cover_set), breakdown)

# ============================================================
# =============  Zsigmondy (p^k - 1) exceptions  =============
# ============================================================

def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n-1)) == 0

def zsigmondy_exception_for_pk_minus_1(p: int, k: int) -> bool:
    if k == 2 and is_power_of_two(p+1):
        return True
    return False

# ============================================================
# =====================  Choix du r  =========================
# ============================================================

def choose_primitive_r_streaming(p: int, e: int, y: int) -> Tuple[Optional[int], Optional[int], str]:
    """Choisit un r | σ(p^e) avec ord_r(p) | (e+1), en factorisant 'juste ce qu'il faut'.
       Stratégie : si un r < y existe => renvoie-le pour L1 immédiate.
       Sinon, parmi r ≥ y, choisit la “balle empoisonnée” (Sophie, grand r, tau(r-1) petit).
    """
    k = e + 1
    if zsigmondy_exception_for_pk_minus_1(p, k):
        return (None, None, "zsigmondy:k=2_and_p+1_is_power_of_two")
    s = sigma_p_pow_e(p, e)
    fac = factorint(s)  # rapide via Pollard–Rho + cache
    cand_lt = []
    cand_ge = []
    for r in fac.keys():
        if r == p:
            continue
        if math.gcd(p, r) != 1:
            continue
        t = order_mod(p, r)
        if t is None:
            continue
        if k % t != 0:
            continue
        if r < y:
            cand_lt.append((r, t))
        else:
            # score “balle empoisonnée”
            sophie = is_probable_prime((r - 1)//2) if (r - 1) % 2 == 0 else False
            tau_r1 = tau_of(r - 1)
            cand_ge.append((r, t, sophie, tau_r1))
    if cand_lt:
        # L1 direct : n'importe lequel suffit (choisissons le plus grand < y)
        r, t = max(cand_lt, key=lambda x: x[0])
        return (r, t, "force_L1")
    if not cand_ge:
        return (None, None, "no_primitive_found")
    # poison : préf. Sophie, puis grand r, puis tau(r-1) petit
    r, t, sophie, tau_r1 = max(cand_ge, key=lambda z: (z[2], z[0], -z[3]))
    return (r, t, "ok_ge_y_poison")

# ============================================================
# ====================  BnB structures  ======================
# ============================================================

@dataclass
class Edge:
    src: int
    dst: int
    order_t: int  # t | (e_src + 1)

@dataclass
class BranchState:
    y: int
    omega_cap: int
    support: Dict[int,int] = field(default_factory=dict)  # p -> e (0 = pas fixé)
    euler_prime: Optional[int] = None                    # unique base à exposant impair
    G_edges: List[Edge] = field(default_factory=list)    # dépendances q -> r
    reason: Optional[str] = None
    v2_sum: int = 0                                      # somme v2(σ(p^e)) pour e>0

    def omega(self) -> int:
        return len(self.support)

    def add_edge(self, q: int, r: int, t: int):
        self.G_edges.append(Edge(q, r, t))

    def has_cycle(self) -> bool:
        g = defaultdict(list)
        nodes = set(self.support.keys())
        for e in self.G_edges:
            g[e.src].append(e.dst)
            nodes.add(e.src); nodes.add(e.dst)
        indeg = defaultdict(int)
        for u in nodes:
            indeg[u] = 0
        for u in g:
            for v in g[u]:
                indeg[v] += 1
        dq = deque([u for u in nodes if indeg[u]==0])
        seen = 0
        while dq:
            u = dq.popleft()
            seen += 1
            for v in g[u]:
                indeg[v] -= 1
                if indeg[v] == 0:
                    dq.append(v)
        return seen != len(nodes)

# ============================================================
# ==============  UB restant (affiné, O(1))  =================
# ============================================================

def S_partial(support: Dict[int,int]) -> float:
    val = 1.0
    for p, e in support.items():
        if e > 0:
            val *= S_factor(p, e)
    return val

def local_X_threshold(state: BranchState) -> int:
    # seuil local = max(y, max(support)) — on peut raffiner si besoin
    if not state.support:
        return state.y
    return max(state.y, max(state.support.keys()))

# ============================================================
# ================   Coupes L0–L3 optimisées   ===============
# ============================================================

def cut_L0_v2(state: BranchState) -> bool:
    # v2(σ(n)) doit être = 1 en final ; si >1, coupe immédiate
    return state.v2_sum > 1

def cut_L1(r: int, y: int) -> bool:
    return r is not None and r < y

def cut_L2_cycle(state: BranchState) -> bool:
    return state.has_cycle()

def neutralite_q_adique_impossible(state: BranchState, UB_room: int) -> Optional[int]:
    incoming: Dict[int,List[Tuple[int,int]]] = defaultdict(list)
    for e in state.G_edges:
        if e.dst in state.support and e.src in state.support and state.support[e.src] > 0:
            incoming[e.dst].append((e.src, state.support[e.src]))
    for q, e_q in state.support.items():
        if e_q == 0 or q == 2:
            continue
        contribs = incoming.get(q, [])
        v_lb = 0
        for (p, e_p) in contribs:
            v_lb += v_q_of_sigma_p_e(q, p, e_p)
        if v_lb < e_q and UB_room == 0:
            return q
    return None

# ============================================================
# ================  Heuristiques (Axe 2)  ====================
# ============================================================

def prime_risk_score(p: int) -> Tuple[int, int]:
    # score = (tau(p-1), -p) — favorise bcp de diviseurs, puis p petit
    return (tau_of(p - 1), -p)

def pick_base_to_expand(state: BranchState) -> Optional[int]:
    # choisir parmi bases avec e = 0 (pas encore fixé)
    candidates = [p for p, e in state.support.items() if e == 0]
    if not candidates:
        # sinon, prendre une base récente (heuristique simple)
        if state.support:
            return min(state.support.keys())
        return None
    # “maillon faible” : plus grand score
    return max(candidates, key=prime_risk_score)

def sort_exponents_for_branch(e_list: List[int]) -> List[int]:
    # prioriser e avec beaucoup de diviseurs pour e+1, puis e grand
    return sorted(e_list, key=lambda e: (tau_of(e + 1), e), reverse=True)

# ============================================================
# =====================  BnB (cœur)  =========================
# ============================================================

@dataclass
class BnBConfig:
    y: int
    omega_cap: int
    even_exps: List[int]
    euler_exps: List[int]
    max_nodes: int = 1000
    sieve_max: int = 200000
    ub: Optional[UBPrimes] = None

@dataclass
class BnBNodeLog:
    node_id: int
    parent_id: Optional[int]
    action: str
    support_snapshot: Dict[int,int]
    cut: Optional[str]
    S_part: float
    UB_rest: float
    S_times_UB: float

def build_ub_helper(y: int, omega_cap: int, sieve_max: int) -> UBPrimes:
    # Primes via crible, puis liste ≥ y (impairs), et préfixe UB
    allp = sieve_primes(sieve_max)
    ge_y = [p for p in allp if p >= y and p % 2 == 1]
    ub = UBPrimes(start_y=y, sieve_max=sieve_max, primes_list=allp, primes_ge_y=[], ub_prefix=[])
    prod = 1.0
    for p in ge_y[:max(1, omega_cap)]:
        if p % 2 == 0:
            continue
        ub.primes_ge_y.append(p)
        prod *= (p/(p-1))
        ub.ub_prefix.append(prod)
    if not ub.primes_ge_y:
        ub.primes_ge_y = []
        ub.ub_prefix = []
    return ub

def UB_rest_affine(cfg: BnBConfig, state: BranchState) -> float:
    used = state.omega()
    m = max(0, cfg.omega_cap - used)
    X = local_X_threshold(state)
    return cfg.ub.ub_from_X(X, m)

def branch_on_base(state: BranchState, cfg: BnBConfig, parent_id: Optional[int], logs: List[BnBNodeLog], next_id: List[int]) -> None:
    if state.omega() > cfg.omega_cap:
        state.reason = "omega_exceeds_cap"
        return

    # L0 : v2 coupe immédiate
    if cut_L0_v2(state):
        logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id, action="L0_v2_cut",
                               support_snapshot=dict(state.support), cut="L0_v2>1",
                               S_part=S_partial(state.support),
                               UB_rest=UB_rest_affine(cfg, state),
                               S_times_UB=S_partial(state.support)*UB_rest_affine(cfg, state)))
        next_id[0] += 1
        return

    # Déficience rapide
    S = S_partial(state.support)
    UB = UB_rest_affine(cfg, state)
    SxUB = S * UB
    logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id, action="check",
                           support_snapshot=dict(state.support),
                           cut=None, S_part=S, UB_rest=UB, S_times_UB=SxUB))
    next_id[0] += 1
    if SxUB < 2.0:
        state.reason = "deficiency_cut"
        return

    # Choisir base à étendre (heuristique)
    p = pick_base_to_expand(state)
    if p is None:
        state.reason = "no_base_to_expand"
        return

    # Exposants autorisés pour p
    allowed_exps = []
    if state.euler_prime is None:
        allowed_exps.extend(cfg.euler_exps)  # impairs ≡ 1 (mod 4)
        allowed_exps.extend(cfg.even_exps)
    elif state.euler_prime == p:
        allowed_exps.extend(cfg.euler_exps)
    else:
        allowed_exps.extend(cfg.even_exps)
    allowed_exps = sort_exponents_for_branch(sorted(set(allowed_exps)))

    # Branches (évalue d'abord les e “explosifs”)
    for e in allowed_exps:
        child = BranchState(y=state.y, omega_cap=state.omega_cap,
                            support=dict(state.support),
                            euler_prime=state.euler_prime,
                            G_edges=list(state.G_edges),
                            v2_sum=state.v2_sum)

        # Marquer Euler si e impair
        if e % 2 == 1:
            if child.euler_prime is None or child.euler_prime == p:
                child.euler_prime = p
            else:
                # Conflit Euler -> L2
                logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                       action=f"assign e={e} to p={p} => bad_euler",
                                       support_snapshot=dict(child.support), cut="L2_euler_conflict",
                                       S_part=S_partial(child.support),
                                       UB_rest=UB_rest_affine(cfg, child),
                                       S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
                next_id[0] += 1
                continue

        # Affecter e et v2_sum
        prev_e = child.support.get(p, 0)
        child.support[p] = e
        # delta v2
        if prev_e > 0:
            child.v2_sum -= v2_sigma_p_e(p, prev_e)
        if e > 0:
            child.v2_sum += v2_sigma_p_e(p, e)

        # L0 : v2
        if cut_L0_v2(child):
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id, action=f"set p={p}^e={e} => L0_v2",
                                   support_snapshot=dict(child.support), cut="L0_v2>1",
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # Déficience rapide après mise à jour
        S = S_partial(child.support)
        UB = UB_rest_affine(cfg, child)
        SxUB = S * UB
        logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                               action=f"after assign p={p}^e={e} check",
                               support_snapshot=dict(child.support), cut=None,
                               S_part=S, UB_rest=UB, S_times_UB=SxUB))
        next_id[0] += 1
        if SxUB < 2.0:
            child.reason = "deficiency_cut"
            continue

        # Choisir un r primitif (factorisation tardive + cache)
        r, t, reason = choose_primitive_r_streaming(p, e, cfg.y)
        if r is None:
            # Pas de primitif “propre” (Zsigmondy/no_primitive) -> coupe conservatrice
            child.reason = f"L2_zsigmondy_or_no_primitive ({reason})"
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                   action=f"set p={p}^e={e}",
                                   support_snapshot=dict(child.support), cut=child.reason,
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # L1 : r < y
        if cut_L1(r, cfg.y):
            child.reason = f"L1_r_lt_y (p={p}, e={e}, r={r})"
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                   action=f"set p={p}^e={e}",
                                   support_snapshot=dict(child.support), cut=child.reason,
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # Ajouter r (si absent) avec e=0
        if r not in child.support:
            if child.omega() + 1 > cfg.omega_cap:
                child.reason = "omega_cap_hit_on_add_r"
                logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                       action=f"add r={r} from p={p}^e={e}",
                                       support_snapshot=dict(child.support), cut=child.reason,
                                       S_part=S_partial(child.support),
                                       UB_rest=UB_rest_affine(cfg, child),
                                       S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
                next_id[0] += 1
                continue
            child.support[r] = 0

        # Arête p -> r
        child.add_edge(p, r, t)

        # L2 : cycle ?
        if cut_L2_cycle(child):
            child.reason = "L2_cycle_detected"
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                   action=f"edge {p}->{r} (t|{e+1})",
                                   support_snapshot=dict(child.support), cut=child.reason,
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # L3 : neutralité q-adique impossible ?
        room = cfg.omega_cap - child.omega()
        q_bad = neutralite_q_adique_impossible(child, room)
        if q_bad is not None:
            child.reason = f"L3_neutralite_impossible_at_q={q_bad}"
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                   action=f"edge {p}->{r} (t|{e+1})",
                                   support_snapshot=dict(child.support), cut=child.reason,
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # Budget de nœuds ?
        if next_id[0] >= cfg.max_nodes:
            child.reason = "node_budget_exhausted"
            logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                                   action=f"edge {p}->{r}",
                                   support_snapshot=dict(child.support), cut=child.reason,
                                   S_part=S_partial(child.support),
                                   UB_rest=UB_rest_affine(cfg, child),
                                   S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
            next_id[0] += 1
            continue

        # Continuer
        logs.append(BnBNodeLog(node_id=next_id[0], parent_id=parent_id,
                               action=f"advance p={p}^e={e}, add r={r}",
                               support_snapshot=dict(child.support), cut=None,
                               S_part=S_partial(child.support),
                               UB_rest=UB_rest_affine(cfg, child),
                               S_times_UB=S_partial(child.support)*UB_rest_affine(cfg, child)))
        nid = next_id[0]
        next_id[0] += 1
        branch_on_base(child, cfg, nid, logs, next_id)

# ============================================================
# ===================  Interface CLI  ========================
# ============================================================

def group_int(n: int) -> str:
    return f"{n:,}".replace(",", " ")

def float_s(x: float, d: int=12) -> str:
    return f"{x:.{d}f}"

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

def cmd_table(args):
    ys = parse_ys(args.ys)
    if args.ymin and args.ymax:
        cand = [n for n in range(args.ymin, args.ymax+1) if n>=3 and n%2==1 and is_probable_prime(n)]
        ys = sorted(set(ys + cand))
    if not ys:
        raise SystemExit("Aucun y fourni (utilise --ys ou --ymin/--ymax).")

    rows = []
    for y in ys:
        M = build_M_y(args.kmax, allow_even_exps=True, allow_one_odd_euler=True)
        Dram, breakdown = D_ram_of_y(y, M)
        omega_cap = args.omega_cap if args.omega_cap else 0
        Cpar = max(0, omega_cap - 1) if omega_cap else None
        row = {
            "y": y,
            "K_max": args.kmax,
            "M_y": M,
            "omega_odd_My": len([p for p in factorint(M) if p != 2]),
            "D_ram": Dram,
            "C_par": Cpar
        }
        rows.append(row)

    header = " y | K_max | omega_odd(M_y) | D_ram(y) | C_par(y;ω_cap) |  M(y) "
    print(header)
    print("-"*len(header))
    for r in rows:
        print(f"{r['y']:3d} | {r['K_max']:5d} | {r['omega_odd_My']:14d} | {r['D_ram']:8d} | {str(r['C_par']):>13} | {r['M_y']}")

    if args.csv_out:
        with open(args.csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["y","K_max","M_y","omega_odd_My","D_ram","C_par"])
            for r in rows:
                w.writerow([r["y"], r["K_max"], r["M_y"], r["omega_odd_My"], r["D_ram"], r["C_par"]])
        print(f"\n[OK] CSV écrit : {args.csv_out}")

def cmd_bnb(args):
    y = args.y
    if not (y and y >= 3 and y % 2 == 1 and is_probable_prime(y)):
        raise SystemExit("Arg --y : fournir un premier impair (ex: 11).")

    even_exps = parse_ys(args.even_exps)
    euler_exps = parse_ys(args.euler_exps)
    even_exps = [e for e in even_exps if e >= 2 and e % 2 == 0]
    euler_exps = [e for e in euler_exps if e >= 1 and (e % 2 == 1) and (e % 4 == 1)]

    if not even_exps:
        even_exps = [2]
    if not euler_exps:
        euler_exps = [1]

    even_exps = sorted(set(even_exps))
    euler_exps = sorted(set(euler_exps))

    cfg = BnBConfig(
        y=y,
        omega_cap=args.omega_cap,
        even_exps=even_exps,
        euler_exps=euler_exps,
        max_nodes=args.max_nodes,
        sieve_max=args.sieve_max
    )

    cfg.ub = build_ub_helper(cfg.y, cfg.omega_cap, cfg.sieve_max)

    # état initial : support = {y : 0}
    state0 = BranchState(y=y, omega_cap=cfg.omega_cap, support={y: 0}, v2_sum=0)

    logs: List[BnBNodeLog] = []
    next_id = [1]
    branch_on_base(state0, cfg, None, logs, next_id)

    # Impression du journal
    print("=== BnB+L run (optimisé) ===")
    print(f"y = {y}  |  omega_cap = {args.omega_cap}  |  even_exps = {even_exps}  |  euler_exps = {euler_exps}  |  max_nodes = {args.max_nodes}  |  sieve_max = {args.sieve_max}\n")
    print(" id | parent | action                           | cut                              |   S_part        |   UB_rest       |    S×UB         | support")
    print("-"*150)
    for i, rec in enumerate(logs, 1):
        pid = rec.parent_id if rec.parent_id is not None else 0
        cut = rec.cut if rec.cut else "-"
        supp = " ".join(f"{p}^{e}" for p,e in sorted(rec.support_snapshot.items()))
        print(f"{rec.node_id:4d} | {pid:6d} | {rec.action:<30} | {cut:<30} | {rec.S_part:14.10f} | {rec.UB_rest:14.10f} | {rec.S_times_UB:14.10f} | {supp}")

    # résumé
    n_def = sum(1 for r in logs if r.cut == "deficiency_cut")
    n_L0 = sum(1 for r in logs if r.cut == "L0_v2>1")
    n_L1 = sum(1 for r in logs if r.cut and r.cut.startswith("L1_"))
    n_L2 = sum(1 for r in logs if r.cut and r.cut.startswith("L2"))
    n_L3 = sum(1 for r in logs if r.cut and r.cut.startswith("L3"))
    n_cap = sum(1 for r in logs if r.cut == "omega_cap_hit_on_add_r")
    n_budget = sum(1 for r in logs if r.cut == "node_budget_exhausted")

    print("\n--- Résumé ---")
    print(f"Nodes enregistrés : {len(logs)}")
    print(f"Déficience (S_part × UB_rest < 2) : {n_def}")
    print(f"Coupes L0(v2): {n_L0}  |  L1 : {n_L1}   |  L2 : {n_L2}   |  L3 : {n_L3}")
    print(f"Atteintes du cap ω : {n_cap}")
    print(f"Budget de nœuds épuisé : {n_budget}")

    if args.csv_out:
        with open(args.csv_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id","parent","action","cut","S_part","UB_rest","S_times_UB","support"])
            for rec in logs:
                pid = rec.parent_id if rec.parent_id is not None else 0
                supp = " ".join(f"{p}^{e}" for p,e in sorted(rec.support_snapshot.items()))
                w.writerow([rec.node_id, pid, rec.action, rec.cut or "", f"{rec.S_part:.12f}", f"{rec.UB_rest:.12f}", f"{rec.S_times_UB:.12f}", supp])
        print(f"\n[OK] CSV écrit : {args.csv_out}")

# ============================================================
# ====================   main CLI parser  ====================
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_table = sub.add_parser("table", help="Tableaux M(y), D_ram(y), C_par(y)")
    ap_table.add_argument("--ys", type=str, help="liste de y (premiers impairs) séparés par des virgules")
    ap_table.add_argument("--ymin", type=int, help="borne inf y (incluse)")
    ap_table.add_argument("--ymax", type=int, help="borne sup y (incluse)")
    ap_table.add_argument("--kmax", type=int, default=40, help="borne K_max pour M(y)")
    ap_table.add_argument("--omega-cap", type=int, help="cap ω (pour afficher C_par = ω_cap-1)")
    ap_table.add_argument("--csv-out", type=str, help="fichier CSV de sortie")
    ap_table.set_defaults(func=cmd_table)

    ap_bnb = sub.add_parser("bnb", help="BnB + verrous (optimisé)")
    ap_bnb.add_argument("--y", type=int, required=True, help="premier impair (ex: 11)")
    ap_bnb.add_argument("--omega-cap", type=int, required=True, help="cap sur le # de bases (ω)")
    ap_bnb.add_argument("--even-exps", type=str, default="2", help="liste d'exposants pairs ≥2 (ex: 2,4,6)")
    ap_bnb.add_argument("--euler-exps", type=str, default="1", help="liste d'exposants Euler (impairs ≡1 mod 4)")
    ap_bnb.add_argument("--max-nodes", type=int, default=5000, help="limite du journal / exploration")
    ap_bnb.add_argument("--sieve-max", type=int, default=200000, help="borne du crible initial")
    ap_bnb.add_argument("--csv-out", type=str, help="CSV du journal BnB")
    ap_bnb.set_defaults(func=cmd_bnb)

    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    # Seed pour Pollard–Rho (non déterministe, c'est ok)
    random.seed(0xC0FFEE)
    main()
