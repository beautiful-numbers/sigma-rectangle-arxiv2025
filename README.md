# Sigma Rectangle arXiv 2025 Companion

Companion code and certificates for  
**The Sigma Rectangle Framework and the Non-Existence of Odd Perfect Numbers** DOI: https://zenodo.org/records/17383511

This repository lets the referee regenerate pivot and threshold constants, verify UB, PNE, and MET certificates, and run a short descent demo above the pivot.  
All routines use exact rational arithmetic and emit auditable outputs that match the appendices.

---

## ✅ What’s included

**Orchestration**  
Entry scripts to run the core checks exactly as used in the paper.  
They provide a small demo, plus a curated set of scenarios for quick verification.

- proof_algorithms.py: public routines wired for the proof workflow  
- run_algo.py: command line entry for the small descent demo  
- run_catalog.py: curated scenarios that cover finite windows and large y

**Budget vs Cost pivot and constants**  
Builds the universal thresholds and the quadratic pivot barrier from explicit prime products.  
Compares minimal complexity against available digit budget and expands until the pivot is crossed.

- budget_vs_cost/bnb_lock_engine.py  
- budget_vs_cost/build_constants_table.py  
- budget_vs_cost/calculate_m_requis.py  
- budget_vs_cost/calculate_d_structure.py  
- budget_vs_cost/calculate_prime_density.py  
- budget_vs_cost/compare_budget_cost.py   
- budget_vs_cost/expand_until_pivot.py  
- budget_vs_cost/E48_table_small.csv  

**Core routines and runs UB PNE MET — src**  
Implements multiplicative UB cuts, p-adic non-equality rules, and the exact MET protocol.  
Includes reproducible runs for ceiling bounds, exact digits control, rough profile reads, and abundance thresholds.

- src/blockA.py  
- src/blockB.py  
- src/blockC.py  
- src/run_ceiling_rigorous.py  
- src/run_exact_digits.py  
- src/run_rough_profile.py  
- src/run_seuil_abondance.py  
- src/A_D6.json  
- src/B_D6.json  
- src/C_D6.json  

**Certificates samples — certs**  
Certificate files that correspond to the formats used in the paper.  
They can be checked directly with the verification commands below.

- certs/certificats_MET_ancres.txt  
- certs/certificats_B4_branchI.txt  
- certs/certificats_B4_auto_anchor.txt
- certs/certificats_B4_omega7.txt 
- certs/appendix_C_supplement.json

**Paper**  
DOI: https://zenodo.org/records/17383511


---
## Quickstart

Requires Python 3.11+. No external packages.

### 1) Effective constants — fix c1 and Y0
Recomputes the effective constants table as in the appendix. Prints ε(y), H(y..y²), RS bound, Φ(y), then derives c1 and Y0.
```bash
python budget_vs_cost/build_constants_table.py --ys 607,613,641,701,997,1499
# write a CSV as well:
python budget_vs_cost/build_constants_table.py --ys 607,613,641,701,997,1499 --csv-out constants_sample.csv
# or compute on a prime range:
python budget_vs_cost/build_constants_table.py --ymin 607 --ymax 2003 --csv-out constants_range.csv
```

### 2) Pivot / budget crossing — X*(D), D_pivot(y), and window table
Shows the calibrated kappa (quadratic growth constant), prints the crossing X*(D), and a small window around it to locate possible abundance at fixed D.
```bash
# basic run
python run_algo.py

# inspect a specific D, print the window rows
python -c "import run_algo as ra; Xs, rows = ra.window_table(100000); print('X*=', Xs); [print(r) for r in rows]"
```

**Get the pivot at a given y and compare budget vs cost.**  
Use `expand_until_pivot.py` to print `D_pivot` (budget) and a progress table of the observed structural cost.
```bash
# zsigmondy-based lower bound (fast structural LB cost model)
python budget_vs_cost/expand_until_pivot.py --y 79 --scenario B --mode zsig-lb \
  --m-requis 912 --max-depth 64 --zsig-per-node 2 --respect-y-rude --print-header

# factor mode (trial division on sigma(p^e); direct but costlier)
python budget_vs_cost/expand_until_pivot.py --y 101 --scenario A --mode factor \
  --m-requis 1331 --max-depth 32 --cap-sigma-digits 200 --respect-y-rude --print-header
```
- **pivot** → `D_pivot = ceil(m * log10(y))` printed when `--m-requis` is provided.  
- **croissance quadratique** → the calibrated **kappa** reported by `run_algo.py`, used in `X*(D)`.  
- **budget** → set via `--m-requis <m>`.  
- **cost** → choose the model with `--mode {zsig-lb,factor}`.

### 3) Parent-minimal supports — the 27 anchors
Generates the 27 parent-minimal supports of size 8, excluding 3, with {5,7} forced. Prints the count and the list.

macOS / Linux / WSL (bash/zsh)
```bash
python - <<'PY'
from pprint import pprint
from proof_algorithms import generate_minimal_supports
S = generate_minimal_supports()
print("count =", len(S))
pprint(S)
PY
```

Windows (PowerShell or CMD)
```powershell
python -c "from pprint import pprint; import sys; sys.path.append('.'); from proof_algorithms import generate_minimal_supports as g; S=g(); print('count =', len(S)); pprint(S)"
```
---

## What the remaining files are for
The other files remain for reproducibility and deeper audit; they are not required for the Quickstart above.

- run_catalog.py: generates the 27 parent supports and builds deterministic MET anchors; prints a table and writes `appendix_C_supplement.json` and `certificats_MET_ancres.txt` at the repo root.
- src/run_ceiling_rigorous.py, src/run_exact_digits.py, src/run_rough_profile.py, src/run_seuil_abondance.py: targeted local verifications (ceiling bounds, exact-digit control, rough profiles, abundance thresholds).
- budget_vs_cost/expand_until_pivot.py: stepwise expansion to the pivot; prints a progress table; does not write `expand_*` artifacts.
- budget_vs_cost/bnb_lock_engine.py: locking/branch engine used in budget-versus-cost analyses.
- budget_vs_cost/calculate_m_requis.py, budget_vs_cost/calculate_d_structure.py, budget_vs_cost/calculate_prime_density.py: helper computations used to build the constants table and cross-check RS bounds.
- budget_vs_cost/E48_table_small.csv: small auxiliary table; only used if consumed by the scripts.
- src/A_D6.json, src/B_D6.json, src/C_D6.json: single sources of truth for D6 parameters, consumed where needed by downstream scripts (not by Block A).

