# TaxoTox — Technical Methods Reference

> **Purpose**: Citable methods reference for the TaxoTox paper.
> Sections are added incrementally as each feature is implemented.
> Last updated: 2026-07-23

---

## 1. Overview

TaxoTox is an R Shiny application for aquatic mixture toxicity assessment.
It calculates Toxic Units (TU) and the Pollution Toxicity Index (PTI) for chemical
mixtures detected in water samples, using species-sensitivity data from the US EPA
ECOTOX Knowledgebase as the primary toxicological denominator.

The application supports three taxonomic groups: fish, algae, and crustacean.
Default behaviour uses the ECOTOX median LC50/EC50 as the TU denominator (Concentration
Addition framework). Advanced methods — HC5-based TU, Benchmark Hazard Index, Independent Action, and CAMA —
are opt-in via an "Advanced Assessment" sidebar panel. When enabled, calculations run for
all selected methods on the same button press. Results are delivered as two separate Excel
files: `taxotox_basic_*.xlsx` (standard CA results, unchanged from default behaviour) and
`taxotox_advanced_*.xlsx` (one sheet per advanced method × taxonomic group). When both files
are produced they are bundled into a single `.zip` download. This separation intentionally
keeps operational outputs clean while making research-grade comparisons available in Excel.

---

## 2. Reference Table (`taxotox_reference.fst`)

### 2.1 Rationale

Previous versions of TaxoTox computed `median(conc_ng_L)` at runtime from the full
ECOTOX dataset (~1–2 million rows). This approach was computationally expensive and
precluded adding alternative denominators without restructuring the app.

The redesigned architecture pre-aggregates all denominators at install time into a
single reference table (`taxotox_reference.fst`), one row per `cas_number × ecotox_group`.
The app performs a simple column lookup at runtime. The raw ECOTOX dataset
(`final_ecotox_data.fst`) is used only during installation and is not deployed to
the production server.

### 2.2 Schema

| Column | Type | Description |
|---|---|---|
| `cas_number` | character | CAS Registry Number |
| `ecotox_group` | character | `fish` / `algae` / `crustacean` |
| `chemical_name` | character | Canonical name from ECOTOX |
| `n_ecotox` | integer | Number of ECOTOX test results (0 = model-only row) |
| `median_lc50_ng_L` | numeric | Median of all LC50/EC50 values (ng/L); current TU denominator |
| `hc5_ssd_ng_L` | numeric | HC5 from log-normal SSD fitted to ECOTOX data (n ≥ 5) |
| `hc5_model_ng_L` | numeric | HC5 from empirical scaling factor (all rows with median LC50) |
| `hc5_method` | character | `"SSD"` / `"Scaled"` / `"Both"` / `NA` |
| `predicted_lc50_ng_L` | numeric | CompTox/OPERA model prediction (ng/L); fish and crustacean only |
| `lc50_source` | character | `"ECOTOX_median"` / `"CompTox_predicted"` / `"Both"` / `NA` |
| `moa_group` | character | Mode of action group from Kramer et al. (2024); `"unknown"` if unclassified (see Section 8) |
| `moa_source` | character | Source of Mode of Action classification: `"Kramer2024"` / `"unknown_in_Kramer2024"` / `"not_in_Kramer2024"` |
| `benchmark_usepa_fish_acute_ng_L` | numeric | US EPA — freshwater fish, acute (ng/L) |
| `benchmark_usepa_crust_acute_ng_L` | numeric | US EPA — freshwater invertebrate, acute (ng/L) |
| `benchmark_usepa_algae_acute_ng_L` | numeric | US EPA — nonvascular plant (algae), acute (ng/L) |
| `benchmark_eu_eqs_aa_marine_ng_L` | numeric | EU Water Framework Directive — Annual Average Environmental Quality Standard for marine/estuarine waters (ng/L) |
| `stc_nowell_ng_L` | numeric | Nowell et al. (2014) sensitive toxicity concentration, read directly from their published Appendix B table (ng/L); crustacean row = cladoceran-genus-restricted only (see Section 10.6) |
| `n_stc_nowell` | integer | Nowell et al.'s published "No. bioassays" count behind `stc_nowell_ng_L` |

---

## 3. Toxicity Data Sources

### 3.1 ECOTOX Knowledgebase

**Source**: US EPA ECOTOX Knowledgebase (https://cfpub.epa.gov/ecotox/).
Downloaded via the `ECOTOXr` R package (`ECOTOXr::download_ecotox_data()`); stored
locally as a SQLite file and updated approximately quarterly.

**Cache versioning and updates**: `download_ecotox_data()` never overwrites or deletes
a previous release — each download adds a new dated
`ecotox_ascii_MM_DD_YYYY.sqlite` (plus a sibling extracted-ASCII folder, `.log`, and
`_cit.txt` file) to the local ECOTOXr cache directory
(`ECOTOXr::get_ecotox_path()`). `taxotox_install.R`'s Step 0a checks for a newer
release: interactively, it prompts before downloading; non-interactively (e.g.
scheduled/automated runs), it checks the release version via
`ECOTOXr::check_ecotox_version()` and downloads automatically only when a newer
release actually exists, so unattended installs stay current without ever
re-downloading data that's already up to date (a failed/offline version check fails
open and just uses the cached database). `get_ecotox_sqlite_file()` — used both by
Step 1 below and by `validate_benchmarks.R` — always auto-selects the most recently
dated release present in the cache, so no code changes are ever needed after an
update. Step Z, the final step of `taxotox_install.R`, cleans up superseded releases
once a rebuild against the newest one has completed: prompted interactively, or
automatically (guarded by a row-count sanity check against the pre-run baseline) on
non-interactive runs.

**Query logic** (`taxotox_install.R`, Steps 1–4):

- Taxa: fish (`ecotox_group LIKE '%fish%'`), algae (`'%algae%'`), crustacean (`'%crustacean%'`)
- Endpoints: LC50, EC50 (and their log-transformed variants `(log)LC50`, `(log)EC50`;
  decorated variants `LC50*`, `LC50/` — decorators stripped after back-transformation)
- Concentration: `min_concentration` = lowest of up to three reported concentration columns
  (`conc1_mean`, `conc2_mean`, `conc3_mean`); conservative choice for mixture toxicity

**Unit normalisation** (Step 2):

All concentrations converted to ng/L using the following factors:

| ECOTOX unit | Alias | Factor to ng/L |
|---|---|---|
| ng/L | ppt | 1 |
| µg/L | ppb | 1 × 10³ |
| mg/L | ppm | 1 × 10⁶ |
| g/L | 0/00 | 1 × 10⁹ |

Rows with other units (e.g. nmol/L, % body weight) are excluded.
The `AI` (active ingredient) prefix is stripped before unit matching.

**Endpoint back-transformation**: entries reported as `(log)LC50` or `(log)EC50`
are back-transformed via `10^value` before unit conversion.

**Aggregation** (Step 7a): median of all `conc_ng_L` values per
`(cas_number, ecotox_group)` pair → `median_lc50_ng_L`. Note that LC50 and EC50
are combined before computing the median; see Section 5.3 for justification.

---

## 4. Standard TU/PTI Method

### 4.1 Formula

$$TU_{ij} = \frac{C_{ij}}{LC50_i}$$

$$PTI_j = \sum_i TU_{ij}$$

where $C_{ij}$ is the measured concentration of compound $i$ in sample $j$ (ng/L),
and $LC50_i$ is the species-sensitivity denominator for compound $i$ from the reference table.

### 4.2 Denominator priority

| Priority | Condition | Denominator used |
|---|---|---|
| 1 | `median_lc50_ng_L` is not NA | ECOTOX median |
| 2 | `predicted_lc50_ng_L` is not NA | CompTox/OPERA prediction |

The denominator is selected at runtime using `dplyr::coalesce(median_lc50_ng_L, predicted_lc50_ng_L)`:
the first non-NA value in priority order is used. Compounds with neither value available
contribute no TU and do not appear in output tables.

Previous versions of TaxoTox computed `median(conc_ng_L)` at runtime by scanning the full
raw ECOTOX dataset (~118,000 rows). Step A-2 replaced this with a direct column lookup from
`taxotox_reference.fst`, where the median is pre-computed at install time. The TU values
produced are identical (confirmed by V-4), and startup time is substantially reduced because
the large raw dataset is no longer loaded.

Compounds using a predicted denominator are flagged with a `~` prefix in output table
column headers and Excel export sheets (e.g. `~ibuprofen`). The flag is applied to any
compound where `median_lc50_ng_L` is NA and `predicted_lc50_ng_L` was used instead,
across all taxonomic groups. Plot axis labels are not flagged to preserve readability.

### 4.3 PTI risk thresholds

| PTI range | Interpretation |
|---|---|
| ≥ 1.0 | Acute risk (mixture concentration exceeds LC50 equivalent) |
| 0.1 – 1.0 | Chronic risk |
| < 0.1 | Low risk |

---

## 5. HC5 Methods

### 5.1 Background and Formula

The 5th percentile of the Species Sensitivity Distribution (HC5, "Hazardous Concentration
for 5% of species") is a regulatory standard for deriving Environmental Quality Standards
(ANZG 2018; EU TGD, 2011). Using HC5 as the TU denominator yields a more protective
assessment than the median LC50, as it targets the most sensitive 5% of species rather
than the median-sensitive species.

$$TU_{ij}^{HC5} = \frac{C_{ij}}{HC5_i}$$

$$PTI_j^{HC5} = \sum_i TU_{ij}^{HC5}$$

where $HC5_i$ is selected by priority: `hc5_ssd_ng_L` if available (fitted SSD),
otherwise `hc5_model_ng_L` (scaling factor estimate). The interpretation thresholds
(PTI ≥ 1.0 acute, 0.1–1.0 chronic) are unchanged from the standard method, but
because HC5 ≪ median LC50, HC5-based PTI values are systematically higher — more
compounds will exceed thresholds.

Two separate HC5 estimates are retained in `taxotox_reference.fst` so the source is
always traceable in output and audit trails:

| Column | Method | When populated |
|---|---|---|
| `hc5_ssd_ng_L` | Log-normal SSD fit to ECOTOX data | n_ecotox ≥ 5 (~15–20% of rows) |
| `hc5_model_ng_L` | Empirical scaling from median LC50 | All rows with any LC50 value |
| `hc5_method` | `"SSD"` / `"Scaled"` / `"Both"` | Always populated if either HC5 exists |

### 5.2 SSD-based HC5 (`hc5_ssd_ng_L`) — step I-1a

A log-normal Species Sensitivity Distribution is fitted to all LC50/EC50 values for
each `(cas_number, ecotox_group)` pair with at least 5 data points, using the
`ssdtools` R package (Thorley & Schwarz, 2018):

```r
ssd_fit_dists(data.frame(Conc = values), dists = "lnorm")
ssd_hc(fit, proportion = 0.05, ci = FALSE)
```

**Minimum n = 5 requirement**: regulatory frameworks (ANZG 2018; EU TGD) require data
from at least 5 species representing at least 4 taxonomic groups for a valid SSD.
Within TaxoTox, each `ecotox_group` (fish, algae, crustacean) is treated as a single
group for SSD fitting, so the 5-data-point threshold serves as a practical lower bound
rather than a full multi-group validation. Rows with fewer than 5 data points receive
NA for `hc5_ssd_ng_L` and fall back to the scaling-factor estimate (Section 5.3).

Fitting failures (non-convergence, degenerate data) return NA silently.

Coverage: approximately 15–20% of compound × group pairs in the current ECOTOX extract
have n_ecotox ≥ 5 and a valid SSD.

### 5.3 Scaling-factor HC5 (`hc5_model_ng_L`) — step I-1b

Standard regulatory practice (EU TGD, REACH) applies a fixed Assessment Factor (AF)
when fewer than 5 species are available — typically AF = 10 for a 3-taxon dataset,
meaning HC5_estimated = median_LC50 / 10.

TaxoTox replaces the fixed AF with an **empirically calibrated scaling factor** derived
from data-rich compounds in ECOTOX:

$$k = \text{median}\left(\frac{\text{hc5\_ssd\_ng\_L}}{\text{median\_lc50\_ng\_L}}\right) \quad \text{per ecotox\_group}$$

$$\text{hc5\_model\_ng\_L} = \text{median\_lc50\_ng\_L} \times k$$

$k$ is computed separately for fish, algae, and crustacean because inter-species
variability — and therefore the HC5/median ratio — differs between taxonomic groups.

The scaling factor is saved to `Data/hc5_scaling_factors.csv` for audit and reuse
without re-running the full installation.

**Comparison to EU TGD AF approach**: both approaches derive a protective estimate from
a single-species or few-species dataset. The key difference is that the EU fixed AF = 10
is a regulatory default based on historical experience, whereas TaxoTox's $k$ is
calibrated to the actual ECOTOX distribution for each taxonomic group. A $k$ consistent
with 0.1 (i.e. AF ≈ 10) would confirm alignment with EU guidance; deviations indicate
that the empirical distribution for a given group warrants a different protection factor.

**Note on CompTox gap-fill rows** (n_ecotox = 0, added at step I-2): for compounds with
no ECOTOX data, `hc5_model_ng_L = predicted_lc50_ng_L × k`. This assignment is deferred
until `predicted_lc50_ng_L` is populated.

### 5.4 Combining LC50 and EC50 before computing the median

TaxoTox merges LC50 (lethal concentration) and EC50 (effective concentration) results
before computing `median_lc50_ng_L`. Both endpoint types represent the concentration
causing a 50% response in an acute toxicity test and are considered interchangeable for
the TU framework in standard practice (Backhaus & Faust, 2012; ECHA guidance).

This assumption is validated empirically by Plot 3 of `Code/validate_reference.R`,
which plots median LC50 vs median EC50 per compound × group for all pairs where both
endpoint types are available. A close 1:1 relationship supports the combined-median
approach; systematic divergence would indicate the need to separate them.

---

## 6. Benchmark Hazard Index

### 6.1 Rationale

The standard Toxic Unit / Pollution Toxicity Index method uses the ECOTOX median
lethal concentration (LC50) as the denominator, which represents the concentration
affecting 50% of the tested species. National regulatory water quality guidelines
use a different protection philosophy: they set benchmarks at concentrations intended
to protect a defined proportion of species (typically 95–99%) under long-term or
acute exposure scenarios.

The Benchmark Hazard Index (also called Hazard Quotient sum) uses these regulatory
benchmarks as denominators instead of LC50 values:

$$HQ_i = \frac{C_i}{BM_i}$$

$$HI = \sum_i HQ_i$$

where $C_i$ is the measured concentration of compound $i$ and $BM_i$ is the selected
national benchmark. The interpretation thresholds are the same as for Pollution Toxicity
Index: $HI \geq 1$ indicates exceedance of the benchmark, $HI < 0.1$ indicates low risk.

The method currently implements two frameworks, US EPA and EU EQS (see Section 6.3 for
why two other national frameworks TaxoTox previously supported were removed, and why EU
EQS specifically was removed and then reinstated).

### 6.2 Benchmark Source

All benchmark values are stored in `taxotox_reference.fst` as columns in ng/L.
Source files are in `Data/` and are parsed by `Code/taxotox_install_benchmarks.R`.
Source values are in µg/L; multiplication by 1000 converts to ng/L.

#### US EPA Aquatic Life Benchmarks

**Source**: US EPA Office of Pesticide Programs — Aquatic Life Benchmarks.
**File**: `Data/USEPA_aquatic_benchmarks.csv`

The US EPA table provides benchmarks for three taxonomic groups at two protection levels
(acute and chronic). The column names in the source file are duplicated; R auto-renames
them with `.1` suffixes. Column mapping:

| R column after read | Taxon | Protection level | Reference table column |
|---|---|---|---|
| `AcuteA` | Freshwater fish | Acute | `benchmark_usepa_fish_acute_ng_L` |
| `AcuteA.1` | Freshwater invertebrate | Acute | `benchmark_usepa_crust_acute_ng_L` |
| `AcuteC` | Nonvascular plant (algae) | Acute | `benchmark_usepa_algae_acute_ng_L` |

**Denominator selection**: acute columns only. The Hazard Quotient framework divides
measured concentrations by an acute protective value; using chronic benchmarks as Hazard
Index denominators would mix protection levels and invalidate the HI interpretation.
Chronic columns are present in the source file but are not used as denominators.

#### EU WFD Environmental Quality Standards

**Source**: EU Water Framework Directive, Directive 2013/39/EU Annex I (~45 priority
substances).
**File**: `Data/EUEPA_aquatic_benchmarks.csv`

Unlike the US EPA table, this source provides only a single Annual Average EQS value per
substance (no separate acute/chronic or taxon-specific split) — `benchmark_eu_eqs_aa_
marine_ng_L`. Where a substance is listed under a parent/sub-name grouping (e.g. aldrin
and dieldrin reported under a combined cyclodiene pesticides entry), the sub-name row's
CAS number is preferred.

### 6.3 Coverage, and Why Two Frameworks Were Removed (and EU EQS Reinstated)

TaxoTox previously also implemented Australia/New Zealand ANZG 2018 Default Guideline
Values at the 95% species protection level (`benchmark_au_anzg_fresh_ng_L` /
`_marine_ng_L`) and Canada CCME freshwater long-term guidelines
(`benchmark_ca_ccme_fresh_lt_ng_L`). Both were removed from
`Code/taxotox_install_benchmarks.R` and the app's benchmark sub-selector because their
compound coverage in TaxoTox's reference table was too low to be useful:

| Framework | Coverage (of matched compound rows) |
|---|---|
| US EPA (fish / crustacean / algae) | ~10.4% / ~10.4% / ~8.6% |
| EU EQS | ~0.9% |
| AU ANZG (freshwater / marine) | ~2.8% / ~1.3% |
| CA CCME | ~2.1% |

EU EQS was removed for the same low-coverage reason initially, but reinstated after
review: its coverage as a *percentage* of TaxoTox's full reference table looks similar to
the other removed frameworks, but that denominator includes TaxoTox's entire, much
broader compound universe (CompTox gap-fill compounds, Nowell-published compounds, etc.).
In absolute terms the EU list's ~45 priority substances were judged a useful enough panel
to keep as a second framework, unlike AU ANZG/CA CCME.

For compounds absent from a selected framework's table, the Hazard Quotient for that
compound is not computed and contributes zero to that framework's Hazard Index — this
means the Hazard Index underestimates risk when coverage is poor for the compounds in a
given sample. Coverage is reported at the end of `Code/taxotox_install_benchmarks.R`'s
console output on every run.

The two removed frameworks' source files (`Data/Australia_aquatic_benchmarks.csv`,
`Data/Canada_aquatic_benchmarks.csv`) were left in `Data/` in case a future revision of
those sources improves coverage enough to reconsider — re-adding them would mean
restoring the I-8/I-9 parsing logic described in `Docs/TaxoTox_Redesign_Plan.md`'s version
history for this file.

### 6.4 References

- US EPA (various years). *Aquatic Life Benchmarks and Ecological Risk Assessments for
  Registered Pesticides*. Office of Pesticide Programs, Washington DC.
- European Union (2013). *Directive 2013/39/EU of the European Parliament and of the
  Council amending Directives 2000/60/EC and 2008/105/EC as regards priority substances
  in the field of water policy*. Annex I: Environmental Quality Standards.

---

## 7. CompTox/OPERA LC50 Gap-Fill

### 7.1 Rationale

ECOTOX covers approximately 8,000–10,000 distinct compounds across the three target
taxonomic groups, but `Known_CAS.fst` — the curated list of compounds monitored in
TaxoTox deployments — includes compounds for which no experimental aquatic toxicity
data exist in ECOTOX (pharmaceuticals, PFAS, industrial chemicals of emerging concern).
Without a denominator, these compounds cannot contribute to PTI or TU calculations.

To maximise coverage without compromising data quality, TaxoTox uses the US EPA
CompTox Chemicals Dashboard (Williams et al., 2017) to retrieve OPERA-model predicted
LC50 values for compounds absent from or data-poor in ECOTOX. OPERA (OPEn structure–
activity/property Relationship App) is an open-source QSAR platform developed and
maintained by the EPA (Mansouri et al., 2018) with demonstrated regulatory acceptance
for chemical hazard screening.

Predicted values are stored in a separate column (`predicted_lc50_ng_L`) and are
never mixed with ECOTOX-derived medians in the same column. The `lc50_source` flag
in every reference row makes the data origin explicit in all outputs.

### 7.2 API

**Service**: EPA CompTox Chemicals Dashboard CTX REST API  
**Base URL**: `https://comptox.epa.gov/ctx-api`  
**Authentication**: free API key via `x-api-key` request header  
**Registration**: https://www.epa.gov/comptox-tools

### 7.3 Workflow (step I-2 in `taxotox_install.R`)

**Step I-2b — CAS → DTXSID resolution**

Each CAS Registry Number in `Known_CAS.fst` is resolved to a DSSTox Substance
Identifier (DTXSID) via a per-compound GET request:

```
GET /chemical/search/equal/{casrn}
```

The first exact CASRN match is used; if none, the top-ranked result is taken.
The DTXSID is the stable identifier used for all subsequent property lookups.
Note: the search endpoint does not return molecular weight.

**Step I-2b2 — Molecular weight (batch)**

Molecular weight (`averageMass`, g/mol) is retrieved via a separate batch call:

```
POST /chemical/detail/search/by-dtxsid/
Body: JSON array of DTXSIDs (max 1000 per request)
```

MW is needed for unit conversion in step I-2d and is not available from the
search or predicted-property endpoints.

**Step I-2c — OPERA predictions (batch)**

Predicted ecotoxicity properties are fetched via:

```
POST /chemical/property/predicted/search/by-dtxsid/
Body: JSON array of DTXSIDs (max 1000 per request)
```

Properties extracted:

| propName | Organism | Endpoint | ecotox_group |
|---|---|---|---|
| `96 Hour Fathead Minnow LC50` | *Pimephales promelas* | LC50 | fish |
| `48 Hour Daphnia Magna LC50` | *Daphnia magna* | LC50 | crustacean |

**Algae**: no suitable OPERA endpoint is available. The only alternative
(*Tetrahymena pyriformis* IGC50) is a ciliate protozoan, not an alga, and is
scientifically inappropriate as a substitute. Predicted LC50 for algae is left as
NA; algal TU can only be computed if ECOTOX data exist.

**Step I-2d — Unit conversion**

The CTX API returns predicted values in mol/L (confirmed empirically; not −log10).
Conversion to ng/L using MW from step I-2b2:

$$C_{\text{ng/L}} = C_{\text{mol/L}} \times MW_{\text{g/mol}} \times 10^9$$

### 7.4 Integration into reference table

- **Existing ECOTOX rows**: `predicted_lc50_ng_L` is added as an additional column.
  `lc50_source` is set to `"Both"` if both values are present, `"ECOTOX_median"` if only
  ECOTOX data exist.
- **New gap-fill rows** (`n_ecotox = 0`): created for Known_CAS compounds with no ECOTOX
  entry but a valid CompTox prediction. `lc50_source = "CompTox_predicted"`.
  `hc5_model_ng_L = predicted_lc50_ng_L × k` (group scaling factor from step I-1b).
- **Algae gap-fill rows**: not created — no predicted denominator available.

### 7.5 References

- Mansouri, K., et al. (2018). OPERA models for predicting physicochemical properties
  and environmental fate endpoints. *Journal of Cheminformatics*, 10(1), 10.
  https://doi.org/10.1186/s13321-018-0263-1

- Williams, A.J., et al. (2017). The CompTox Chemistry Dashboard: a community data
  resource for environmental chemistry. *Journal of Cheminformatics*, 9(1), 61.
  https://doi.org/10.1186/s13321-017-0247-6

---

## 8. Mode of Action Classification

### 8.1 Role in the CAMA Method

The Concentration Addition with Mode of Action grouping (CAMA) method partitions
compounds in a sample into mechanistic groups, applies Concentration Addition within
each group, and then combines the group-level effect estimates using Independent Action
across groups (see Section 9). The partition relies on a **Mode of Action group** assigned
to every compound in the reference table (`moa_group` column).

In standard Concentration Addition (used for the Toxic Unit / Pollution Toxicity Index
calculation), all compounds are assumed to act by the same non-specific mechanism — an
assumption known to be violated for mixtures of compounds with distinct modes of action.
CAMA relaxes this assumption while remaining computationally tractable.

### 8.2 Mode of Action Groups

`moa_group` values come directly from an external curated database (Section 8.3) rather
than a fixed, hardcoded list — the actual set of groups present in `taxotox_reference.fst`
depends on which compounds are in that source's coverage. A representative sample of
`moa_group` values (verified directly against the source data during development):

| `moa_group` value | Mechanism | Example compounds |
|---|---|---|
| `Neuromuscular system` | Includes acetylcholinesterase (AChE) inhibition — accumulation of acetylcholine at nerve synapses | Chlorpyrifos, chlorpyrifos-methyl, diazinon (organophosphates) |
| `Photosynthesis inhibition` | Blocks photosynthetic electron transport (largely Photosystem II) in plants and algae | Atrazine, diuron (triazine and phenylurea herbicides) |
| `Neuroactive` | Includes sodium channel modulation and other neural mechanisms | Permethrin and other pyrethroids |
| `unknown` | No Mode of Action classification available | Any compound not covered, or not yet classified, by the source database |

Some entries are themselves multi-mechanism, comma-joined labels straight from the source
data (e.g. `"Cardiovascular system, Neuromuscular system"`) — these are kept verbatim
rather than split, since CAMA only needs a consistent group label to partition compounds
for Concentration Addition, not a single "primary" mechanism.

Compounds with no classifiable information are assigned `moa_group = "unknown"` — an
explicit, honest label, not a default assumption of narcosis/baseline toxicity. See
Section 8.4 for why this changed from an earlier narcosis-default approach, and Section
10.5 for how CAMA handles the `"unknown"` group in its output.

### 8.3 Classification Strategy

Mode of Action group assignment is implemented in `Code/taxotox_install_moa.R`. It runs
automatically at the end of `taxotox_install.R` (Step I-12, alongside
`taxotox_install_benchmarks.R` and `taxotox_install_nowell.R` — see Section 10.6 for the
latter), but is also fully runnable standalone afterward, e.g. to refresh just the MoA
join without re-running the rest of the pipeline. Source: Kramer, L., Schulze, T.,
Klüver, N., Altenburger, R.,
Hackermüller, J., Krauss, M., & Busch, W. (2024). Curated mode-of-action data and effect
concentrations for chemicals relevant for the aquatic environment. *Scientific Data*, 11,
26. https://doi.org/10.1038/s41597-023-02904-7 — a systematically curated, peer-reviewed,
freely downloadable dataset (Zenodo, DOI 10.5281/zenodo.10071824, no login required) from
the Helmholtz Centre for Environmental Research (UFZ), covering fish, algae, and
crustaceans separately — the same three taxonomic groups TaxoTox itself uses.

The dataset ships as two CSV tables (`Data/Kramer2024_TableB_chemicals.csv`,
`Data/Kramer2024_TableC_moa.csv`), both keyed by an internal `ID` and joined together on
that key. Table B supplies `cas_rn` (standard hyphenated CAS, converted to TaxoTox's
dash-free `cas_number` convention); Table C supplies `MoA_broad`, used verbatim as
`moa_group`. Compounds with no CAS number in the source (`cas_rn == "n/a"`, mostly
unregistered mixtures/generic entries) are dropped before joining, since they cannot be
matched to `taxotox_reference.fst` by CAS. The join was confirmed 1:1:1 (no duplicate
`ID` in either table, no duplicate `cas_rn` in Table B) during development, so no
tie-breaking logic is needed.

`moa_source` records provenance:

- `"Kramer2024"` — the compound has a real MoA classification in the source data
- `"unknown_in_Kramer2024"` — the compound is in the source dataset, but its own
  `MoA_broad` value is `"unknown"` there (the source's own way of marking a chemical it
  could not classify)
- `"not_in_Kramer2024"` — the compound isn't in the source dataset at all (it covers
  ~3,387 chemicals relevant to the aquatic environment, not TaxoTox's full compound
  universe, which also includes CompTox-gap-filled and Nowell-published compounds)

Both `"unknown_in_Kramer2024"` and `"not_in_Kramer2024"` result in `moa_group = "unknown"`
— TaxoTox does not distinguish "the source couldn't classify it" from "the source doesn't
have it" in the group label itself, only in `moa_source`.

### 8.4 Coverage and Why the Narcosis Default Was Removed

An earlier version of this classification used a ~100-name regex pattern list for known
pesticide classes (organophosphates, carbamates, triazines, phenylureas, pyrethroids),
falling back to a LogKow-based Verhaar narcosis check that contributed no classifications
in practice, and finally defaulting **every unmatched compound to `"narcosis"`** as a
conservative placeholder. In the dataset that heuristic ran against, this meant ~97% of
compound × group rows were narcosis by default — not because 97% of compounds were
actually confirmed to act by baseline toxicity, but because the classifier had no better
information for them.

The Kramer et al. (2024) dataset was adopted specifically because it distinguishes real
narcosis/baseline-toxicity classifications from compounds it genuinely could not classify
(labelled `"unknown"` in the source itself). Re-running the join against TaxoTox's current
reference table gives, per unique compound:

| `moa_source` | Compounds |
|---|---|
| `Kramer2024` (real classification) | 1,702 |
| `unknown_in_Kramer2024` | 594 |
| `not_in_Kramer2024` | 4,096 |

Real MoA coverage (~27% of TaxoTox's ~6,400 unique compounds) is far lower than the old
scheme's ~97% narcosis figure — by design. The old figure measured how often the
classifier assigned *some* label; this one measures how often that label reflects an
actual curated classification.

Key limitations of the current classification:

1. **Coverage is bounded by the source dataset's scope**: Kramer et al. (2024) covers
   ~3,387 chemicals "relevant for the aquatic environment" — a broad but not exhaustive
   set. Compounds outside that scope (including many CompTox-gap-filled and
   Nowell-published compounds already in `taxotox_reference.fst`) get `moa_group =
   "unknown"` regardless of whether they have a real, documented mode of action elsewhere
   in the literature.

2. **`"unknown"` is not a risk-neutral default**: unlike the old narcosis-default
   approach (deliberately conservative, since narcosis is a weak/common mechanism),
   `"unknown"` compounds form their own CAMA group and are combined with all other groups
   via Independent Action. Section 10.5 covers how this affects `E_mix`, and why CAMA now
   produces a second output variant that excludes the unknown-MoA group entirely.

3. **Mode of Action versus mode of toxic action**: as before, the classification used
   here is the toxicological mode of action at the organism level (how the chemical kills
   or inhibits), not the molecular initiating event in the sense of the Adverse Outcome
   Pathway (AOP) framework.

### 8.5 References

- Kramer, L., Schulze, T., Klüver, N., Altenburger, R., Hackermüller, J., Krauss, M., &
  Busch, W. (2024). Curated mode-of-action data and effect concentrations for chemicals
  relevant for the aquatic environment. *Scientific Data*, 11, 26.
  https://doi.org/10.1038/s41597-023-02904-7 — dataset:
  https://doi.org/10.5281/zenodo.10071824

- Verhaar, H.J.M., van Leeuwen, C.J., Hermens, J.L.M. (1992). Classifying environmental
  pollutants. *Chemosphere*, 25(4), 471–491.
  https://doi.org/10.1016/0045-6535(92)90280-5

---

## 9. Validation of the Reference Table

### 9.1 Overview

The reference table is validated by `Code/validate_reference.R`, which produces an
interactive HTML report (`Docs/validate_reference.html`). The report covers four areas:

1. HC5 computation pipeline (Plots 1–4)
2. LC50/EC50 endpoint interchangeability (Plot 3)
3. CompTox/OPERA gap-fill calibration (Plot 5, when applicable)
4. Mode of Action classification coverage (Plot 7)
5. Benchmark framework coverage per compound (Plot 8)

A targeted spot-check script (`Code/temp_v2_spotcheck.R`) prints reference table values
for five well-studied compounds against published reference ranges. It can be run at any
time after installation and deleted after review.

### 9.2 HC5 Spot-Check Compounds (V-2)

Five compounds were selected for manual validation of `median_lc50_ng_L` and
`hc5_ssd_ng_L` against published values. All values in ng/L.

| Compound | CAS | Groups | Expected median LC50 (ng/L) | Expected HC5 (ng/L) | Reference |
|---|---|---|---|---|---|
| Atrazine | 1912-24-9 | fish, algae | 100,000 – 10,000,000 | 1,000 – 20,000 | ECOTOX; ANZG (2018) |
| Chlorpyrifos | 2921-88-2 | fish, crustacean | 10 – 1,000 | 0.1 – 1 | ECOTOX; ANZG (2018) marine |
| Diuron | 330-54-1 | fish, algae | 1,000,000 – 100,000,000 | 500 – 5,000 | ECOTOX; EU EQS basis |
| Permethrin | 52645-53-1 | fish, crustacean | 100 – 10,000 | 10 – 100 | ECOTOX; literature |
| Copper | 7440-50-8 | fish, crustacean | 10,000 – 1,000,000 | 100 – 2,000 | ECOTOX; ANZG (2018) marine |

**Flag criterion**: a value more than 10× outside the expected range is flagged
`<<< CHECK` by the spot-check script. Values within 10× are considered acceptable
given the natural variability of toxicity data across laboratories and species.

**Note on CAS format**: ECOTOX stores CAS numbers without hyphens (e.g. `"191224-9"`
is stored as `"1912249"`). The spot-check script strips hyphens before lookup.

### 9.3 Published Reference Sources for Spot-Check

- **ANZG (2018)**: Australian and New Zealand Guidelines for Fresh and Marine Water
  Quality. Table 3.4.1 (freshwater) and Table 3.4.2 (marine) — Default Guideline Values
  at 95% species protection (LOSP 95). Available at:
  https://www.waterquality.gov.au/anz-guidelines/guideline-values/default/water-quality-toxicants

- **EU EQS basis documents**: European Commission Technical Guidance for Deriving
  Environmental Quality Standards (TGD-EQS, 2011). Substance-specific basis documents
  available via the EU Common Implementation Strategy (CIS) website:
  https://circabc.europa.eu/ui/group/8ae6f840-4a26-4c82-a26a-3614a16b4d5d

- **ECOTOX Knowledgebase**: US EPA ECOTOX Knowledgebase (https://cfpub.epa.gov/ecotox/).
  Published summary statistics for median LC50 per compound × taxonomic group are
  consistent with ranges reported in Sanderson et al. (2003) and Posthuma et al. (2002).

- **Permethrin literature**: Laskowski (2002) — *A dose-response synthesis for selected
  pyrethroid insecticides*. Acute LC50 for rainbow trout: 2–5 µg/L; Daphnia magna: 0.6–2 µg/L.
  Converted: fish LC50 ~ 2,000–5,000 ng/L; crustacean ~ 600–2,000 ng/L.

### 9.4 Validation Status

| Step | Description | Status |
|---|---|---|
| V-1 | Run `taxotox_install.R` end-to-end; confirm `taxotox_reference.fst` produced | **Confirmed** |
| V-2 | Spot-check HC5 and median LC50 for 5 well-studied compounds | **Confirmed** |
| V-3 | Verify CompTox predictions against published OPERA values | **Not applicable** — see below |
| V-4 | Load `sample_GB.xlsx`; confirm standard TU/PTI results unchanged vs. deployed version | **Confirmed** |

**V-1 result**: `taxotox_reference.fst` produced with 8,365 rows (3 taxonomic groups ×
~2,788 CAS numbers). All planned columns present including Mode of Action and all four
benchmark frameworks.

**V-2 result**: Spot-check confirmed. For all five compounds, `median_lc50_ng_L` and
`hc5_ssd_ng_L` fell within the expected order-of-magnitude range for at least one
taxonomic group. Apparent out-of-range values in remaining groups are expected and
correct: national water quality guidelines (ANZG, EU Environmental Quality Standards)
are derived from the most sensitive taxonomic group in the dataset. The fish and
crustacean rows for algae-sensitive compounds (e.g. atrazine, diuron) and the fish rows
for crustacean-sensitive compounds (e.g. chlorpyrifos, permethrin) will always sit
substantially above the guideline value — this reflects genuine inter-group differences
in sensitivity, not a computation error. The `<<< CHECK` flag in `temp_v2_spotcheck.R`
fires at 10× outside range; values in less-sensitive groups routinely exceed this
threshold by design.

**V-4 result**: The refactored app (`app.R` with `taxotox_reference.fst` denominator lookup)
was run against `Data/sample_GB.xlsx` and its output compared to the previously deployed
version of TaxoTox (which computed `median(conc_ng_L)` at runtime from `final_ecotox_data.fst`).
Pollution Toxicity Index values were identical across all samples and all three taxonomic
groups. The only observable difference was the order of compound columns in the toxicity
worksheets of the Excel export, which changed because the reference table has a different
internal row ordering than the raw ECOTOX dataset. Column order does not affect any
computed value and is not scientifically meaningful.

**V-3 result**: Not applicable. The CompTox gap-fill step (I-2) adds rows only for
compounds with no ECOTOX data (`n_ecotox = 0`). There are therefore no compound × group
pairs with both a measured ECOTOX median and a CompTox prediction, making a direct
empirical comparison impossible with this dataset. The predictions are accepted on the
basis of published OPERA model performance (Mansouri et al., 2018: typical RMSE ≈ 0.6–0.8
log-units on external test sets) and are flagged in all outputs via the `lc50_source` column.

---

## 10. Independent Action (IA)

### 10.1 Background

Independent Action (IA), also called Response Addition or the Hewlett–Plackett model,
is an alternative mixture toxicity framework to Concentration Addition. While CA assumes
all compounds act by the same mechanism and their effects are fully additive at the dose
level, IA assumes compounds act independently — each causes an effect according to its
own dose-response curve, and the probability of the mixture causing no effect is the
product of the probabilities of each individual compound causing no effect.

IA is more appropriate than CA when compounds act via dissimilar mechanisms, particularly
in environmentally realistic mixtures containing both specific-receptor-acting compounds
and baseline narcotics. For the standard TaxoTox monitoring dataset (predominantly narcosis
compounds), CA and IA produce similar results. Divergence increases when the mixture is
dominated by one or two potent specific-acting compounds (e.g. a pyrethroid and an
organophosphate at high TU fractions).

### 10.2 Formula

TaxoTox implements IA with a fixed Hill slope of n = 1 (log-logistic dose-response),
which corresponds to a simple Langmuir saturation curve:

$$E_i = \frac{C_i}{LC50_i + C_i} \tag{4}$$

LaTeX source:
```latex
E_i = \frac{C_i}{LC50_i + C_i}
```

where $E_i$ is the fractional effect (0–1) of compound $i$ at measured concentration $C_i$,
and $LC50_i$ is the same denominator used in the standard TU calculation
(`dplyr::coalesce(median_lc50_ng_L, predicted_lc50_ng_L)`).

The mixture fractional effect is then:

$$E_{mix} = 1 - \prod_{i=1}^{n} (1 - E_i) \tag{5}$$

LaTeX source:
```latex
E_{mix} = 1 - \prod_{i=1}^{n} (1 - E_i)
```

$E_{mix}$ is dimensionless and bounded 0–1, where $E_{mix} = 0.5$ corresponds to a
mixture effect equivalent to the LC50 — the same threshold at which $PTI = 1$ in
the CA framework.

### 10.3 Output

Results are written to three sheets in the advanced Excel file (`IA Algae`, `IA Crustacean`,
`IA Fish`). Each sheet has:

- **`E_mix`** (first column) — mixture fractional effect per sample (0–1)
- Per-compound **`E_i`** columns — individual fractional effects; useful for identifying
  the dominant contributors to mixture risk

### 10.4 Assumptions and Limitations

- **Hill slope n = 1**: the sigmoidal steepness parameter is fixed at 1 for all compounds.
  This is a common simplification in environmental mixture assessment (Backhaus & Faust, 2012).
  Compounds with steeper dose-response curves (n > 1, e.g. some pyrethroids) will have their
  individual fractional effects underestimated at low concentrations and overestimated at high.

- **Missing concentrations treated as zero**: compounds present in the reference table but
  absent from a given sample contribute $E_i = 0$ to the product, leaving $E_{mix}$ unchanged.

- **Denominator same as standard TU**: IA does not require a separate denominator column.
  This means compounds flagged with `~` (CompTox gap-fill) carry the same uncertainty
  into IA calculations as into the standard PTI.

---

## 10.5 CAMA — Concentration Addition with Mode of Action Grouping

### 10.5.1 Background

CAMA (Concentration Addition within Mode of Action groups, Independent Action between groups)
is a hybrid method that combines CA and IA to address the limitations of applying either
framework alone to a multi-mechanistic mixture (Drescher & Boedeker, 1995; Altenburger et al., 2000).

The reasoning: within a group of compounds sharing the same mode of action, CA is the
correct model (joint similarity assumption holds). Between groups with different modes of
action, IA is more appropriate (mechanisms are independent). CAMA exploits this by
applying CA within mechanistic groups first, then IA across groups.

### 10.5.2 Formula

**Step 1 — Within-group CA** (per MoA group $g$, per sample $j$):

$$group\_TU_{gj} = \sum_{i \in g} \frac{C_{ij}}{LC50_i}$$

**Step 2 — Group fractional effect** (logistic transformation):

$$E_{group,gj} = \frac{group\_TU_{gj}}{1 + group\_TU_{gj}}$$

This converts the unbounded CA sum (TU) into a fractional effect (0–1) using the same
Hill n = 1 assumption as IA. A group TU of 1.0 → $E_{group} = 0.5$.

**Step 3 — Between-group IA**:

$$E_{mix,j} = 1 - \prod_g (1 - E_{group,gj})$$

### 10.5.3 Mode of Action Groups

Groups come from an external curated database (see Section 8 for the full classification
methodology) rather than a fixed list — `moa_group` is Kramer et al. (2024)'s own
`MoA_broad` value for each compound (e.g. `Neuromuscular system`, `Photosynthesis
inhibition`, `Neuroactive`), or the literal group `unknown` for compounds with no MoA
classification available (from any cause — see Section 8.3 for the distinction between
`moa_source` values, all of which collapse to `moa_group = "unknown"`).

### 10.5.4 Output

Results are written to **six** sheets — two variants (unknown-MoA compounds included vs.
excluded) for each of the three taxonomic groups:

- `CAMA Algae`, `CAMA Crustacean`, `CAMA Fish` — all matched compounds included, with
  `unknown`-MoA compounds forming their own group like any other
- `CAMA Algae (known MoA)`, `CAMA Crustacean (known MoA)`, `CAMA Fish (known MoA)` —
  compounds with `moa_group == "unknown"` excluded entirely before the within-group CA /
  between-group IA calculation. (Named `"(known MoA)"` rather than the more verbose
  `"(MoA-known only)"` because Excel worksheet names are capped at 31 characters --
  `"CAMA Crustacean (MoA-known only)"` is 32 and openxlsx rejects it.)

Each sheet has:

- **`E_mix`** (first column) — mixture fractional effect per sample (0–1)
- Per-group **`E_group`** columns (e.g. `E_Neuromuscular system`,
  `` `E_Photosynthesis inhibition` ``, `E_unknown`) — useful for seeing which mechanistic
  group drives the mixture risk

Comparing the two variants for the same sample shows how much the unknown-MoA group is
driving the result: since IA treats groups as independent risks that combine
(`E_mix = 1 - ∏(1 - E_group)`), a group with a large `E_group` — including an `unknown`
group dominated by a handful of high-concentration but unclassified compounds — can
noticeably raise `E_mix` even though it says nothing about a shared mechanism. The
"known MoA" variant strips that out, showing what CAMA would report using only
compounds with an actual curated Mode of Action.

### 10.5.5 Practical expectation

Given that real MoA classifications currently cover roughly a quarter of TaxoTox's
compound universe (Section 8.4), most compounds fall into the `unknown` group. This means
the "all compounds" CAMA variant will often behave similarly to standard CA-based PTI for
water samples dominated by unclassified compounds — not because those compounds are known
to act non-specifically (the old narcosis-default assumption), but because CAMA has no
mechanistic information to partition them further. The "known MoA" variant is the
more mechanistically meaningful one whenever a sample contains enough classified
compounds to populate multiple real MoA groups; comparing it against the "all compounds"
variant and against standard PTI is the recommended way to interpret CAMA output.

### 10.5.6 References

- Drescher, K., & Boedeker, W. (1995). Assessment of the combined effects of substances:
  the relationship between concentration addition and independent action. *Biometrics*, 51, 716–730.

- Altenburger, R., et al. (2000). Predictability of the toxicity of multiple chemical
  mixtures to Vibrio fischeri: mixtures composed of similarly acting chemicals.
  *Environmental Toxicology and Chemistry*, 19(10), 2341–2347.

---

## 10.6 Taxon-Sensitive PTI (Nowell et al. 2014)

### 10.6.1 Background

This method replicates the USGS NAWQA Pesticide Toxicity Index (PTI) as defined by
Nowell, Norman, Moran, Martin & Stone (2014), the methodology behind the published
Fish-, Cladoceran-, and Benthic Invertebrate-PTI values reported across USGS NAWQA/RSQA
water-quality studies (e.g. Covert et al., 2020, NWQN water years 2013–2017). It exists
in TaxoTox specifically so results are directly comparable to that published literature
— unlike the Standard PTI (median ECOTOX LC50/EC50, pooled across a whole taxonomic
group) or the Benchmark Hazard Index (a cross-invertebrate EPA regulatory value), both
of which use denominators computed over a different, broader species population than
Nowell et al.'s taxon-restricted approach.

TaxoTox uses Nowell et al.'s own **published** STC values directly (from their
Supplementary Appendix B), rather than re-deriving them from TaxoTox's local ECOTOX
snapshot. An earlier version of TaxoTox recomputed STC from ECOTOX using Nowell et al.'s
"Approach B" procedure (5th percentile of individual toxicity values when n > 12,
otherwise the minimum), but this diverged substantially from the published numbers —
e.g. Diazinon fish STC: ~45 µg/L recomputed vs. 85 µg/L published; Atrazine fish STC:
~2099 µg/L recomputed vs. 4500 µg/L published — because Nowell et al.'s STC also draws
on USEPA OPP aquatic-life benchmarks and the Pesticide Properties Database (PPDB), not
ECOTOX alone, and because ECOTOX itself has been revised since their 2012–2013 data
pull. Using the published table instead makes TaxoTox's Taxon-Sensitive PTI exactly
reproduce, rather than approximate, the values reported in downstream USGS literature
(e.g. Covert et al. 2020).

### 10.6.2 Formula

$$TU_i = \frac{C_i}{STC_i} \qquad PTI_j = \sum_i TU_{ij}$$

Standard Concentration Addition, identical in form to the Standard PTI — only the
denominator source differs. $STC_i$ is read directly from Nowell et al.'s published
Table B.1 (fish) / Table B.2 (cladocerans); see Nowell et al. (2014) Supplementary
Appendix A, sections 3–4, for how they derived it (5th percentile of individual
toxicity values when n > 12, otherwise the minimum — "Approach B" of five candidate
procedures they evaluated).

### 10.6.3 Taxon definition and data source

Source file: `Data/Nowell2014_AppB.xlsx`, the official supplementary Appendix B from
Nowell et al. (2014), downloaded from
https://water.usgs.gov/nawqa/pnsp/pubs/Nowell2014_STOTEN_PTI/Nowell2014_SuppInfo_PTI.zip
and joined into `taxotox_reference.fst` by `Code/taxotox_install_nowell.R` (step I-10).

- **Fish**: `Table B.1 - Fish` (~480 compounds after dropping rows with no assigned CAS
  number).
- **Cladoceran** (reported under the "crustacean" taxonomic group for compatibility with
  TaxoTox's existing per-group calculation code): `Table B.2 - Cladocerans` (~463
  compounds), already restricted by Nowell et al. to the 17 cladoceran genera they
  tested (their Appendix C "Table C.1 PTI taxa") — not TaxoTox's full crustacean group
  (which spans amphipods, copepods, mysids, and other non-cladoceran crustaceans).
- **Algae**: not available. Nowell et al.'s method covers fish, cladocerans, and benthic
  invertebrates only — it was never defined for algae.
- **Benthic invertebrates**: `Table B.3 - Benthic invertebrates` exists in the source
  file but is not joined into TaxoTox — it spans aquatic insects (e.g. midges, mayflies)
  that TaxoTox's ECOTOX extraction does not currently pull in at all (only
  fish/algae/crustacean `ecotox_group`s are queried; see Section 3), and TaxoTox does
  not currently expose a Benthic Invertebrate PTI method in the UI.

**"Non-standard" flagged rows are used as published, not excluded.** Nowell et al.'s
table flags some rows' toxicity source ("Toxicity value type/source" column, e.g.
"non-std OPP", "non-std PPDB") when the underlying test didn't meet their standard
duration/endpoint criteria. TaxoTox does not filter these out: hundreds of compounds
have no standard-duration ECOTOX data at all, so their published STC *is* a non-std
value — excluding it would delete the compound from coverage entirely, not make the
remaining data more standard. The flag documents a known bias direction (Nowell et al.
note whether it over- or under-estimates toxicity); it is not an instruction to discard
the row, and downstream users of the table (e.g. Covert et al. 2020) apply it as-is.

Compounds present in Nowell et al.'s table that TaxoTox's own ECOTOX pull never
surfaced (because their STC came from OPP/PPDB rather than ECOTOX, or because ECOTOX
coverage has shifted since 2012–2013) are added to `taxotox_reference.fst` as new rows,
the same way step I-2 (CompTox gap-fill) adds rows for compounds absent from ECOTOX.

### 10.6.4 Output

Two sheets: `Nowell PTI (Fish)`, `Nowell PTI (Cladoceran)` — column `PTI` (matching
Nowell/Covert's own terminology, distinct from the Benchmark HI method's `HI` column),
plus per-compound TU columns.

### 10.6.5 References

- Nowell, L.H., Norman, J.E., Moran, P.W., Martin, J.D., & Stone, W.W. (2014). Pesticide
  Toxicity Index — a tool for assessing potential toxicity of pesticide mixtures to
  freshwater aquatic organisms. *Science of the Total Environment*, 476–477, 144–157.
  https://doi.org/10.1016/j.scitotenv.2013.12.088
- Covert, S.A., Shoda, M.E., Stackpoole, S.M., & Stone, W.W. (2020). Pesticide mixtures
  show potential toxicity to aquatic life in U.S. streams, water years 2013–2017.
  *Science of the Total Environment*, 745, 141285.
  https://doi.org/10.1016/j.scitotenv.2020.141285

---

## 11. CASRN Resolution

Compound names in the input file are resolved to CAS Registry Numbers (CASRNs) through two sequential layers:

**Layer 1 — Known_CAS (exact, instant):** a curated internal table (`Known_CAS.fst`) of compound names and their CASRNs. Exact string matches are confirmed automatically with no user interaction. The table is designed to hold multiple synonyms per compound (e.g. "Albuterol", "Salbutamol", and "Albuterol / Salbutamol" may all map to the same CASRN 18559-94-9). Each entry has a `PREFERRED_NAME` which is the canonical display name used in all outputs.

**Layer 2 — PubChem REST API (optional):** when enabled, unresolved names are queried live against the PubChem compound database via the `webchem` R package. Candidate matches are presented for user review in the CASRN Matching tab. Accepted matches are appended to a pending-compounds Google Sheet for future integration into Known_CAS via the curation script (`Code/taxotox_curate.R`). Curation is Sheets-only — there is no local-file fallback — so pending items persist across shinyapps.io restarts and there is exactly one queue to review, regardless of which machine or deployment logged the compound.

**Manual entry:** compounds still unresolved after both layers can have their CASRN entered directly in the CASRN Matching tab. These are also logged to the same Google Sheet.

### Synonym handling and duplicate CASRN detection

Known_CAS supports many-to-one name→CASRN mapping to accommodate synonyms and common name variants. If two different rows in the input file resolve to the same CASRN (e.g. the user supplied both "Albuterol" and "Salbutamol" as separate compounds), the app detects this and shows a warning in the CASRN Matching tab. The two rows are processed independently, which means their concentrations are counted separately and summed in the PTI — effectively double-counting the compound. The recommended corrective action is to remove one synonym from the input file before calculating.

---

## 12. Limitations and Caveats

- **Combined LC50/EC50 median**: endpoints are combined on the assumption of
  interchangeability for acute aquatic toxicity. This is assessed empirically
  via `validate_reference.R` Plot 3 but not formally tested per compound.

- **HC5 from scaling factor**: `hc5_model_ng_L` is an order-of-magnitude estimate,
  not a statistically fitted SSD. It should be interpreted as indicative rather than
  regulatory-grade for individual compounds. The group-level $k$ may not reflect the
  inter-species variability of any specific compound.

- **Benchmark coverage for emerging contaminants**: all four national benchmark
  frameworks (US EPA, EU Environmental Quality Standards, Australia ANZG, Canada CCME)
  were developed primarily for pesticides, metals, and legacy pollutants. Coverage of
  pharmaceuticals, personal care products, PFAS, and current-use industrial chemicals —
  which constitute the majority of compounds in a typical estuarine monitoring dataset —
  is sparse to absent. The Hazard Index therefore substantially underestimates risk from
  these compound classes, and absence of a benchmark should not be interpreted as absence
  of hazard. Plot 8 of `validate_reference.R` shows which monitored compounds are covered.

- **CompTox gap-fill predictions are fish and crustacean only**: no suitable OPERA
  endpoint exists for algae. Compounds absent from ECOTOX will have no algal TU
  denominator and will contribute zero to the algal Pollution Toxicity Index regardless
  of their actual algal toxicity.

---

## References

- Backhaus, T., & Faust, M. (2012). Predictive environmental risk assessment of chemical
  mixtures: a conceptual framework. *Environmental Science & Technology*, 46(5), 2564–2573.

- European Chemicals Agency (ECHA). *Guidance on Information Requirements and Chemical
  Safety Assessment, Chapter R.10: Characterisation of dose [concentration]-response for
  environment*. ECHA, Helsinki.

- EU Technical Guidance Document (TGD) on Risk Assessment. Part II. European Commission
  Joint Research Centre, 2003.

- Thorley, J., & Schwarz, C. (2018). ssdtools: An R package to fit Species Sensitivity
  Distributions. *Journal of Open Source Software*, 3(31), 1082.

- ANZG (2018). *Australian and New Zealand Guidelines for Fresh and Marine Water Quality*.
  Australian and New Zealand Governments, Canberra.
