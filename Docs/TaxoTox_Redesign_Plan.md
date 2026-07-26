# TaxoTox Redesign — Multi-Method Toxicity Assessment

**Guiding principle**: default behaviour stays identical for basic users.
All new methods are opt-in under an "Advanced" panel.

---

## New Capabilities

1. **Gap-fill LC50 via CompTox API** — predicted LC50 for compounds absent from ECOTOX
2. **SSD / HC5 denominator** — alternative to ECOTOX median, computed at install time
3. **Benchmark-based Hazard Index** — US EPA and EU EQS (AU ANZG, CA CCME were later
   removed for low coverage; EU EQS was also removed then reinstated — see I-7/I-8/I-9 below)
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
- ~~**Australia ANZG**: LOSP 95% (95% species protection) for both freshwater and marine.~~
- ~~**Canada CCME**: Freshwater Long-Term (chronic) guideline — most comprehensive table.~~
- **Australia ANZG / Canada CCME — removed**: both denominators above were removed from
  `taxotox_install_benchmarks.R` and the app; compound coverage in TaxoTox's reference
  table was too low to be useful (~1.3–2.8%, ~2.1% respectively, vs. US EPA's ~8.6–10.4%).
  See `TaxoTox_Technical_Methods.md` Section 6.3.
- **EU EQS — removed, then reinstated**: initially removed alongside AU ANZG/CA CCME for
  the same low-coverage reasoning (~0.9% of reference-table rows), but restored after
  review: its coverage as an absolute panel (~45 priority substances) was judged useful
  despite the low percentage against TaxoTox's much broader full compound universe. See
  `TaxoTox_Technical_Methods.md` Section 6.3.

---

## TODO

Items are struck through when complete. Updated after each implementation chunk.

### Phase 1 — Install pipeline

- [x] ~~**I-1a** Add SSD/HC5 computation (`ssdtools`, n≥5) to `taxotox_install.R`; populate `hc5_ssd_ng_L`~~
- [x] ~~**I-1b** Compute per-group empirical scaling factor `k = median(hc5_ssd / median_lc50)` from data-rich compounds; apply `hc5_model_ng_L = median_lc50 × k` to all ECOTOX rows; save `k` to `Data/hc5_scaling_factors.csv`~~
- [x] ~~**I-2**  Add CompTox API query for predicted LC50 to `taxotox_install.R`; populate `predicted_lc50_ng_L`~~
- [x] ~~**I-3**  Add CompTox MoA lookup to `taxotox_install.R`; populate `moa_group`, `moa_source`~~ — **superseded**: the name-heuristic + CompTox-LogKow/Verhaar approach defaulted ~97% of compounds to `"narcosis"` with no real classification. Replaced by `Code/taxotox_install_moa.R`, which joins Kramer et al. (2024)'s externally curated Mode-of-Action database instead (real classification for ~1,700 compounds; explicit `"unknown"` for the rest, not a narcosis default). See `TaxoTox_Technical_Methods.md` Section 8.
- [x] ~~**I-4**  Write `taxotox_reference.fst` at end of install (all columns populated)~~
- [x] ~~**I-5**  Write `taxotox_install_benchmarks.R` skeleton~~
- [x] ~~**I-6**  Parse `USEPA_aquatic_benchmarks.csv` — handle duplicate column names; extract acute columns per taxon~~
- [x] ~~**I-7**  Parse `EUEPA_aquatic_benchmarks.csv` — skip 6-row header; extract AA-EQS marine column; strip value qualifiers~~ — **removed, then reinstated**: initially removed for ~0.9% compound coverage in the reference table, but restored to `taxotox_install_benchmarks.R` after review judged its absolute coverage (~45 priority substances) useful despite the low percentage. See `TaxoTox_Technical_Methods.md` Section 6.3.
- [x] ~~**I-8**  Parse `Australia_aquatic_benchmarks.csv` — split Freshwater / Marine rows; use LOSP 95~~ — **removed**: ~1.3–2.8% coverage; parsing logic deleted from `taxotox_install_benchmarks.R`, source CSV left in `Data/` unused.
- [x] ~~**I-9**  Parse `Canada_aquatic_benchmarks.csv` — use Freshwater Long Term; coerce `"No data"` / `"NRG"` → NA~~ — **removed**: ~2.1% coverage; same treatment as I-8.
- [x] ~~**I-10** China GB 3838-2002~~ — excluded (insufficient coverage of emerging contaminants)
- [x] ~~**I-11** Join all parsed benchmarks into `taxotox_reference.fst` by `cas_number`~~ — joins the US EPA (I-6) and EU EQS (I-7) tables after I-8/I-9's removal; the join step itself is no longer separately numbered in `taxotox_install_benchmarks.R`.
- [x] ~~**V-1**  Run `taxotox_install.R` end-to-end; confirm `taxotox_reference.fst` produced correctly~~ 8,365 rows; all planned columns present
- [x] ~~**V-2**  Spot-check `hc5_ssd_ng_L` against published SSDs for 3–5 well-studied compounds~~ Confirmed for atrazine, chlorpyrifos, diuron, permethrin, copper; cross-group CHECK flags expected (guidelines set by most-sensitive group)
- [x] ~~**V-3**  Verify `predicted_lc50_ng_L` (CompTox) matches published OPERA values for 5 test compounds~~ Not applicable — gap-fill targets compounds absent from ECOTOX; no overlap for direct comparison; accepted on basis of published OPERA performance (RMSE ~0.6–0.8 log-units)

### Phase 2 — App core refactor

- [x] ~~**A-1** Replace `final_ecotox_data` startup load with `taxotox_reference`~~
- [x] ~~**A-2** Refactor TU calculation block to use reference table lookup instead of runtime median~~
- [x] ~~**A-3** Add gap-fill badge (`~`) next to compound names where predicted LC50 was used~~
- [x] ~~**V-4** Load `sample_GB.xlsx`; confirm standard TU/PTI results unchanged vs. current version~~ PTI values identical to deployed version; compound column order in toxicity sheets differs (expected — reference table has different internal ordering than raw ECOTOX; values unaffected)

### Phase 3 — Curation script

- [x] ~~**C-1** Write `taxotox_curate.R` interactive loop (read `temp_CAS`, prompt Add/Skip/Reject, write `Known_CAS`)~~
- [x] ~~**C-2** CompTox gap-fill for newly confirmed compounds~~ — **redesigned**: gap-fill is handled by `taxotox_install.R` on the next upgrade run; removed from curate script to avoid duplication. Curate script now only manages Known_CAS; install script handles reference table.

- [x] ~~**V-6** Run `temp_v6_curate_validation.R` to validate `taxotox_curate.R` end-to-end~~ Confirmed: compound removed from Known_CAS, curated via local FST mode, reference rows matched before/after
- [x] **C-3** Curation consolidated to Sheets-only, wired into `taxotox_install.R` — removed the local `temp_CAS.fst` fallback from `taxotox_curate.R` and `app.R` (Google Sheets is now the single source of truth, matching how the app already defaulted `TAXOTOX_SHEET_ID` in practice). Deleted the standalone `update_ecotox.R` entry point (its only content, the ECOTOX-release download, is now Step 0a of `taxotox_install.R`). Added Step 0b: checks the curation Sheet for pending items and offers to run `taxotox_curate.R` — but, like Step 0a, only ever prompts when run interactively; a non-interactive (`Rscript`) run just prints the pending count and continues, since curation needs human judgement per compound and there's no TTY to prompt against during a silent/background install.

### Phase 4 — Advanced features

- [x] ~~**A-4** Add collapsible "Advanced Assessment" panel in sidebar (hidden by default)~~
- [x] ~~**A-5** Add method checkboxes: HC5, Benchmark HI, IA, CAMA~~
- [x] ~~**A-6** Add benchmark sub-selector: US EPA, EU EQS, AU ANZG, CA CCME~~ (appears below all method checkboxes when Benchmark HI is ticked) — **AU ANZG / CA CCME checkboxes later removed** for low coverage; **EU EQS checkbox also removed then reinstated**; sub-selector now shows US EPA and EU EQS.
- [x] ~~**A-7** Advanced Results tab~~ — **design changed**: no screen tab; advanced results go to a second Excel file only
- [x] ~~**A-8** Add coverage warning when selected method has poor compound coverage~~ (amber box in sidebar, threshold < 80%)
- [x] ~~**A-8b** Duplicate CASRN warning~~ — amber box at top of CASRN Matching tab when two input names resolve to the same CASRN
- [x] ~~**A-9** HC5-based TU/PTI calculation~~ — `v$hc5_results`: three sheets (`HC5 Algae/Crustacean/Fish`) written to advanced Excel
- [x] ~~**A-13a** Download logic~~ — single `.xlsx` when basic only; `.zip` with `taxotox_basic_*.xlsx` + `taxotox_advanced_*.xlsx` when advanced methods selected; button labels update accordingly
- [x] ~~**A-10** Benchmark Hazard Index — calculation~~ — `v$benchmark_results`: US EPA produces three group-specific sheets (`Bench. US EPA (Fish/Crust/Algae)`) using acute columns; EU EQS produces one additional sheet (`Bench. EU EQS`); HI first column. (AU ANZG / CA CCME originally each produced one additional sheet too; that code path was removed for low coverage and not reinstated — see I-8/I-9 above.)
- [x] ~~**A-11** IA (Hill n=1) — calculation~~ — `v$ia_results`: three sheets (`IA Algae/Crustacean/Fish`); E_mix first column, per-compound E_i columns after; same LC50 denominator as standard TU
- [x] ~~**A-12** CAMA — MoA grouping → CA → IA~~ — `v$cama_results`: **six** sheets (`CAMA Algae/Crustacean/Fish` + `CAMA Algae/Crustacean/Fish (known MoA)`); E_mix first, then per-MoA E_group columns. Groups come from Kramer et al. (2024) rather than a fixed list (see I-3 above); the "known MoA" variant excludes the `unknown` group entirely, since real classifications currently cover only a fraction of compounds. (Sheet names use `"(known MoA)"` rather than `"(MoA-known only)"` -- the latter pushes `"CAMA Crustacean (MoA-known only)"` to 32 characters, over Excel's 31-character worksheet name limit; this actually broke the advanced download once and was fixed.)
- [x] ~~**A-14** HC5 Method Coverage sheet~~ — **dropped**: coverage is shown in the sidebar warning and in `validate_reference.R` Plot 8; a dedicated Excel sheet adds complexity without new information for the end user
- [x] ~~**A-15** MoA Coverage sheet~~ — **dropped**: MoA coverage is reported in `validate_reference.R` Plot 7; operational users do not need per-compound MoA in the download
- [x] ~~**A-16** Add `lc50_source` column to standard output sheets~~ — **dropped**: gap-fill compounds are already flagged with `~` prefix in column headers; a redundant column would clutter standard output for basic users
- [ ] **V-5** Enable all advanced methods; confirm ZIP download; inspect all sheets

### Phase 5 — Deployment

- [ ] **V-7** Deploy to shinyapps.io; confirm `taxotox_reference.fst` loads and app runs correctly
