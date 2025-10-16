# Sigma Rectangle — arXiv 2025 Companion

Companion code and certificates for  
**“The Sigma Rectangle Framework and the Non-Existence of Odd Perfect Numbers”** (arXiv:25xx.xxxxx).

This repository is a **minimal, referee-oriented** pack to **evaluate the proof**:
- regenerate the **pivot/threshold** constants (Appendix E),
- check **UB** cuts and **PNE/MET** certificates (Section 6),
- run a **small descent** demo above the pivot (Section 8),
- inspect a **sample** of human-readable certificates (Appendix C/H).

> The paper’s TeX source is included as `The_Sigma_Rectangle_Framework_and_the_Non_existence_of_Odd_Perfect_Numbers.tex`.

---

## ✅ What’s included (curated for evaluation)

**Orchestration / API**
- `proof_algorithms.py` — public API for proof routines  
- `run_algo.py` — CLI entrypoint (small descent demo)  
- `run_catalog.py` — curated scenarios (finite window + large-y)

**Budget vs Cost (pivot & constants)**
- `budget_vs_cost/bnb_lock_engine.py`  
- `budget_vs_cost/build_constants_table.py`  
- `budget_vs_cost/calculate_m_requis.py`  
- `budget_vs_cost/calculate_d_structure.py`  
- `budget_vs_cost/calculate_prime_density.py`  
- `budget_vs_cost/compare_budget_cost.py`  *(keep this twin; remove the other)*  
- `budget_vs_cost/expand_until_pivot.py`  
- `budget_vs_cost/E48_table_small.csv`  *(small auxiliary table)*

**Core routines & runs (UB / PNE / MET) — `src/`**
- `src/blockA.py`, `src/blockB.py`, `src/blockC.py`  
- `src/run_ceiling_rigorous.py`  
- `src/run_exact_digits.py`  
- `src/run_rough_profile.py`  
- `src/run_seuil_abondance.py`  
- `src/A_D6.json`, `src/B_D6.json`, `src/C_D6.json`  *(single source of truth for D6 JSONs)*

**Certificates (samples) — `certs/`**
- `certs/certificats_MET_ancres.txt`  
- `certs/certificats_B4_branchI.txt`  
- `certs/certificats_B4_auto_anchor.txt` *or* `certs/certificats_B4_omega7.txt` *(pick one)*  
- `certs/appendix_C_supplement.json`

**Paper**
- `The_Sigma_Rectangle_Framework_and_the_Non_existence_of_Odd_Perfect_Numbers.tex`

---

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
