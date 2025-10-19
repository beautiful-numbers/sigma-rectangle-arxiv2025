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

- proof_algorithms.py — public routines wired for the proof workflow  
- run_algo.py — command line entry for the small descent demo  
- run_catalog.py — curated scenarios that cover finite windows and large y

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
LaTeX source of the manuscript for completeness.  
Reviewers may also use the Zenodo record above for the archived version.


---

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
