# TaxoTox Redesign — Multi-Method Toxicity Assessment

**Guiding principle**: default behaviour stays identical for basic users.
All new methods are opt-in under an "Advanced" panel.

---

## New Capabilities

1. **Gap-fill LC50 via CompTox API** — predicted LC50 for compounds absent from ECOTOX
2. **SSD / HC5 denominator** — alternative to ECOTOX median, computed at install time
3. **Benchmark-based Hazard Index** — US EPA, EU EQS, AU ANZG, CA CCME
4. **Independent Action (IA)** — fixed Hill slope n=1
5. **CAMA** — Mode-of-Action grouped CA within groups → IA between groups
6. **Curation script** — private tool to promote `temp_CAS` → `Known_CAS` and fill reference table

---

## Architecture Note

A new pre-aggregated **`taxotox_reference.fst`** table (one row per `cas_number × ecotox_group`)
replaces the runtime median computation in `app.R`. The app does a simple column lookup.
`final_ecotox_data.fst` (raw test results) stays local — **not deployed to shinyapps.io**.

---

## Method Formulas

| Method | Formula | Denominator |
|---|---|---|
| Standard TU/PTI | TU = C / LC50_median; PTI = Σ TU | `median_lc50_ng_L` (else `predicted_lc50_ng_L`) |
| HC5-based TU/PTI | TU = C / HC5 | `hc5_ssd_ng_L` if available, else `hc5_model_ng_L` |
| Benchmark HI | HQ = C / benchmark; HI = Σ HQ | selected `benchmark_*_ng_L` column |
| IA (Hill n=1) | E_i = C/(LC50+C); E_mix = 1−∏(1−E_i) | `median_lc50_ng_L` |
| CAMA | group_TU = Σ(CA within MoA); E_group = group_TU/(1+group_TU); E_mix = 1−∏(1−E_group) | `median_lc50_ng_L` + `moa_group` |

### HC5 two-column design

Two separate HC5 columns are retained in `taxotox_reference.fst` so source is always traceable:

| Column | Method | When populated |
|---|---|---|
| `hc5_ssd_ng_L` | ECOTOX SSD (log-normal, `ssdtools`) | n_ecotox ≥ 5 (~15–20% of rows) |
| `hc5_model_ng_L` | Empirical scaling: `median_lc50 × k` | all rows with any LC50 value |
| `hc5_method` | "SSD" / "Scaled" / "Both" | always |

The app uses `hc5_ssd_ng_L` preferentially; falls back to `hc5_model_ng_L`.
Both values and the method column are reported in Excel output.

**HC5 scaling factor rationale**: Fitting a log-normal SSD to only 3 CompTox-predicted values
(one per taxonomic group) is not standard practice — regulatory frameworks require ≥5 species
from ≥4 groups (ANZG 2018, EU TGD). Instead, I-1b uses an empirically calibrated scaling factor:

```
k = median(hc5_ssd / median_lc50)   # computed per ecotox_group from n≥5 compounds
hc5_model_ng_L = median_lc50_ng_L × k
```

This is equivalent to the EU TGD Assessment Factor approach (AF ≈ 10 for 3-species data),
but calibrated to the actual ECOTOX distribution rather than a fixed regulatory default.
The scaling factor `k` is saved to `Data/hc5_scaling_factors.csv` for audit/reuse.
For n_ecotox = 0 rows (added at step I-2+), `predicted_lc50_ng_L × k` is used instead.

### Key design decisions documented in code

- **US EPA benchmark denominator**: acute columns only (`AcuteA` fish, `AcuteA.1` invertebrate,
  `AcuteC` algae). The TU formula divides by an acute lethal endpoint (LC50/EC50); using chronic
  benchmarks would mix protection levels. No guideline recommends chronic benchmarks as TU
  denominators within the standard CA framework.
- **EU EQS denominator**: Annual Average EQS for "Other surface waters" (marine/coastal/transitional),
  reflecting TaxoTox's primary use in estuarine and coastal monitoring.
- **Australia ANZG**: LOSP 95% (95% species protection) for both freshwater and marine.
- **Canada CCME**: Freshwater Long-Term (chronic) guideline — most comprehensive table.

---

## TODO

Items are struck through when complete. Updated after each implementation chunk.

### Phase 1 — Install pipeline

- [x] ~~**I-1a** Add SSD/HC5 computation (`ssdtools`, n≥5) to `taxotox_install.R`; populate `hc5_ssd_ng_L`~~
- [x] ~~**I-1b** Compute per-group empirical scaling factor `k = median(hc5_ssd / median_lc50)` from data-rich compounds; apply `hc5_model_ng_L = median_lc50 × k` to all ECOTOX rows; save `k` to `Data/hc5_scaling_factors.csv`~~
- [x] ~~**I-2**  Add CompTox API query for predicted LC50 to `taxotox_install.R`; populate `predicted_lc50_ng_L`~~
- [x] ~~**I-3**  Add CompTox MoA lookup to `taxotox_install.R`; populate `moa_group`, `moa_source`~~
- [x] ~~**I-4**  Write `taxotox_reference.fst` at end of install (all columns populated)~~
- [x] ~~**I-5**  Write `taxotox_install_benchmarks.R` skeleton~~
- [x] ~~**I-6**  Parse `USEPA_aquatic_benchmarks.csv` — handle duplicate column names; extract acute columns per taxon~~
- [x] ~~**I-7**  Parse `EUEPA_aquatic_benchmarks.csv` — skip 6-row header; extract AA-EQS marine column; strip value qualifiers~~
- [x] ~~**I-8**  Parse `Australia_aquatic_benchmarks.csv` — split Freshwater / Marine rows; use LOSP 95~~
- [x] ~~**I-9**  Parse `Canada_aquatic_benchmarks.csv` — use Freshwater Long Term; coerce `"No data"` / `"NRG"` → NA~~
- [x] ~~**I-10** China GB 3838-2002~~ — excluded (insufficient coverage of emerging contaminants)
- [x] ~~**I-11** Join all parsed benchmarks into `taxotox_reference.fst` by `cas_number`~~
- [ ] **V-1**  Run `taxotox_install.R` end-to-end; confirm `taxotox_reference.fst` produced correctly
- [ ] **V-2**  Spot-check `hc5_ssd_ng_L` against published SSDs for 3–5 well-studied compounds
- [ ] **V-3**  Verify `predicted_lc50_ng_L` (CompTox) matches published OPERA values for 5 test compounds

### Phase 2 — App core refactor

- [ ] **A-1** Replace `final_ecotox_data` startup load with `taxotox_reference`
- [ ] **A-2** Refactor TU calculation block to use reference table lookup instead of runtime median
- [ ] **A-3** Add gap-fill badge (`~`) next to compound names where predicted LC50 was used
- [ ] **V-4** Load `sample_GB.xlsx`; confirm standard TU/PTI results unchanged vs. current version

### Phase 3 — Curation script (private — not published; can run in parallel with Phase 2)

- [ ] **C-1** Write `taxotox_curate.R` interactive loop (read `temp_CAS`, prompt Add/Skip/Reject, write `Known_CAS`)
- [ ] **C-2** Add CompTox gap-fill step for newly confirmed compounds
- [ ] **C-3** Add `.gitignore` entry to exclude `taxotox_curate.R` from published repo
- [ ] **V-6** Run `taxotox_curate.R` against a test `temp_CAS.fst`; confirm `Known_CAS` updated correctly

### Phase 4 — Advanced features

- [ ] **A-4** Add collapsible "Advanced Assessment" panel in sidebar (hidden by default)
- [ ] **A-5** Add method checkboxes: HC5, Benchmark HI, IA, CAMA
- [ ] **A-6** Add benchmark sub-selector: US EPA, EU EQS, AU ANZG, CA CCME
- [ ] **A-7** Show "Advanced Results" tab only when ≥1 advanced method is selected
- [ ] **A-8** Add coverage warning when selected method has poor compound coverage
- [ ] **A-9**  HC5-based TU/PTI — calculation + table rendering (prefer `hc5_ssd_ng_L`, fallback `hc5_model_ng_L`)
- [ ] **A-10** Benchmark Hazard Index — calculation + table rendering
- [ ] **A-11** IA (Hill n=1) — calculation + table rendering
- [ ] **A-12** CAMA — MoA grouping → CA → IA + table rendering
- [ ] **A-13** Add advanced method sheets to Excel export:
  - `HC5 Algae`, `HC5 Fish`, `HC5 Crustacean` — include both `hc5_ssd_ng_L` and `hc5_model_ng_L` columns, plus `hc5_method` so source is always visible
  - `Benchmark HI`
  - `IA Results`
  - `CAMA Results`
- [ ] **A-14** Add `HC5 Method Coverage` sheet — table of all compounds showing `n_ecotox`, `hc5_ssd_ng_L`, `hc5_model_ng_L`, `hc5_method`; colour-coded (green = SSD, amber = model, red = neither)
- [ ] **A-15** Add `MoA Coverage` sheet to Excel export
- [ ] **A-16** In standard output sheets, add `lc50_source` column (ECOTOX / CompTox-predicted) so users can see which denominators were empirical vs modelled
- [ ] **V-5** Enable all advanced methods; confirm no crash; inspect IA and CAMA outputs

### Phase 5 — Deployment

- [ ] **V-7** Deploy to shinyapps.io; confirm `taxotox_reference.fst` loads and app runs correctly
