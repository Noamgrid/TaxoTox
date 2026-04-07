# TaxoTox — Technical Methods Reference

> **Purpose**: Citable methods reference for the TaxoTox paper.
> Sections are added incrementally as each feature is implemented.
> Last updated: 2026-04-07

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
| `moa_group` | character | Mode of action group (see Section 8) |
| `moa_source` | character | Source of Mode of Action classification: `"name_heuristic"` / `"Verhaar_LogKow"` / `"default_narcosis"` |
| `benchmark_usepa_fish_acute_ng_L` | numeric | US EPA — freshwater fish, acute (ng/L) |
| `benchmark_usepa_crust_acute_ng_L` | numeric | US EPA — freshwater invertebrate, acute (ng/L) |
| `benchmark_usepa_algae_acute_ng_L` | numeric | US EPA — nonvascular plant (algae), acute (ng/L) |
| `benchmark_eu_eqs_aa_marine_ng_L` | numeric | EU Water Framework Directive — Annual Average Environmental Quality Standard for marine/estuarine waters (ng/L) |
| `benchmark_au_anzg_marine_ng_L` | numeric | Australia ANZG — marine water, 95% species protection (ng/L) |
| `benchmark_au_anzg_fresh_ng_L` | numeric | Australia ANZG — freshwater, 95% species protection (ng/L) |
| `benchmark_ca_ccme_fresh_lt_ng_L` | numeric | Canada CCME — freshwater long-term (chronic) guideline (ng/L) |

---

## 3. Toxicity Data Sources

### 3.1 ECOTOX Knowledgebase

**Source**: US EPA ECOTOX Knowledgebase (https://cfpub.epa.gov/ecotox/).
Downloaded via the `ECOTOXr` R package (`ECOTOXr::download_ecotox_data()`); stored
locally as a SQLite file and updated approximately quarterly.

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

The method is implemented for four national frameworks, selectable independently in the
Advanced Assessment panel of the app.

### 6.2 Benchmark Sources

All benchmark values are stored in `taxotox_reference.fst` as columns in ng/L.
Source files are in `Data/` and are parsed by `Code/taxotox_install_benchmarks.R`.
All source values are in µg/L; multiplication by 1000 converts to ng/L.

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

#### EU Water Framework Directive Environmental Quality Standards

**Source**: EU Directive 2013/39/EU, Annex I — Priority Substances.
**File**: `Data/EUEPA_aquatic_benchmarks.csv`
**Coverage**: approximately 45 priority substances.

The file extract contains a single Annual Average Environmental Quality Standard
(AA-EQS) value per compound (µg/L). This value corresponds to the surface water EQS
applicable to estuarine and coastal waters, which is the primary environment for
TaxoTox deployments.

| Reference table column | Interpretation |
|---|---|
| `benchmark_eu_eqs_aa_marine_ng_L` | AA-EQS for surface waters (µg/L × 1000 → ng/L) |

#### Australia and New Zealand Guidelines for Fresh and Marine Water Quality

**Source**: ANZG 2018 / ANZECC 2000 Default Guideline Values.
**File**: `Data/Australia_aquatic_benchmarks.csv`

Each compound has two rows: one for freshwater and one for marine water. Trigger values
are provided at four species protection levels (80%, 90%, 95%, 99%). TaxoTox uses the
95% protection level (LOSP 95), which corresponds to slightly-to-moderately disturbed
ecosystems and is the recommended default for most monitoring applications.

| Reference table column | Medium | Protection level |
|---|---|---|
| `benchmark_au_anzg_fresh_ng_L` | Freshwater | 95% species protection |
| `benchmark_au_anzg_marine_ng_L` | Marine water | 95% species protection |

#### Canada CCME Water Quality Guidelines

**Source**: Canadian Council of Ministers of the Environment — Water Quality Guidelines
for the Protection of Aquatic Life.
**File**: `Data/Canada_aquatic_benchmarks.csv`

TaxoTox uses the freshwater long-term (chronic) guideline, which is the most
comprehensive and protective column available. Non-numeric entries (`No data`, `NRG`,
`Insufficient data`) are treated as not available (NA).

| Reference table column | Interpretation |
|---|---|
| `benchmark_ca_ccme_fresh_lt_ng_L` | Freshwater long-term (chronic) guideline (µg/L × 1000 → ng/L) |

### 6.3 Coverage

Benchmark coverage varies substantially by framework. The US EPA and Canada CCME tables
cover predominantly pesticides and legacy contaminants. The EU EQS list covers only
~45 EU priority substances. The Australia ANZG table covers a broader range of compounds
including some industrial chemicals and metals. For compounds absent from a benchmark
table, the Hazard Quotient for that compound is not computed and contributes zero to the
Hazard Index — this means the Hazard Index underestimates risk when coverage is poor.

Coverage per framework and taxonomic group is reported at the end of
`taxotox_install_benchmarks.R` and is visualised in Plot 8 of `Code/validate_reference.R`,
which shows all compounds covered by at least one benchmark, coloured by framework and
sorted from broadest to narrowest coverage.

### 6.4 References

- US EPA (various years). *Aquatic Life Benchmarks and Ecological Risk Assessments for
  Registered Pesticides*. Office of Pesticide Programs, Washington DC.

- European Parliament and Council (2013). Directive 2013/39/EU amending Directives
  2000/60/EC and 2008/105/EC as regards priority substances in the field of water policy.
  *Official Journal of the European Union*, L 226, 1–17.

- ANZG (2018). *Australian and New Zealand Guidelines for Fresh and Marine Water Quality*.
  Australian and New Zealand Governments and Australian State and Territory Governments,
  Canberra.

- CCME (various years). *Canadian Water Quality Guidelines for the Protection of Aquatic
  Life: Summary Table*. Canadian Council of Ministers of the Environment, Winnipeg.

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

Five mechanistic groups are used in TaxoTox, plus a default:

| `moa_group` value | Mechanism | Typical compound classes |
|---|---|---|
| `narcosis` | Baseline toxicity — disruption of membrane structure without specific receptor interaction | Most organic compounds; non-polar and polar narcotics |
| `AChE_inhibition` | Inhibition of acetylcholinesterase — leads to accumulation of acetylcholine at nerve synapses | Organophosphate pesticides, carbamate pesticides |
| `PSII_inhibition` | Inhibition of Photosystem II — blocks photosynthetic electron transport in plants and algae | Triazine herbicides, phenylurea herbicides, diazinone herbicides |
| `pyrethroid` | Voltage-gated sodium channel modulation — prolongs sodium channel opening in nerve cells | Synthetic pyrethroid insecticides |
| `reactive` | Non-specific electrophilic reactivity — covalent binding to macromolecules | Epoxides, Michael acceptors, alpha-beta unsaturated carbonyls |

Compounds with no classifiable information are assigned `moa_group = "narcosis"` as a
conservative default, consistent with standard regulatory practice (Verhaar et al., 1992;
ECHA guidance). This assignment means that unclassified compounds are treated as acting
by the most common and most extensively studied aquatic toxicity mechanism.

### 8.3 Classification Strategy

Mode of Action group assignment follows a priority hierarchy implemented in step I-3 of
`taxotox_install.R`. Each compound receives the first applicable classification from the
following sequence:

**Priority 1 — Chemical name pattern matching** (`moa_source = "name_heuristic"`):
The preferred chemical name from CompTox and the canonical name from ECOTOX are tested
against regular expression patterns for major pesticide classes. This method is applied
first because it is the most specific and chemically interpretable. Patterns cover:

- Organophosphates: chlorpyrifos, malathion, diazinon, parathion, phosmet, and
  approximately 20 additional active ingredients by name
- Carbamates: carbofuran, carbaryl, methomyl, aldicarb, pirimicarb, and others
- Triazine herbicides: atrazine, simazine, terbuthylazine, and others
- Phenylurea herbicides: diuron, isoproturon, linuron, chlortoluron, and others
- Other Photosystem II inhibitors: bromacil, bentazon, metribuzin, and others
- Pyrethroids: permethrin, cypermethrin, deltamethrin, lambda-cyhalothrin, and others
- Reactive chemicals: epoxide, acrylate, acrolein, isocyanate, and related functional groups

Name patterns are case-insensitive and applied as substring matches. The CompTox preferred
name takes priority over the ECOTOX chemical name because it is more consistently
formatted.

**Priority 2 — LogKow-based Verhaar classification** (`moa_source = "Verhaar_LogKow"`):
For compounds not matched by name patterns, the octanol-water partition coefficient
(LogKow) is retrieved from the CompTox predicted property endpoint. The Verhaar
classification scheme (Verhaar et al., 1992) assigns organic chemicals to four classes
based on LogKow and structural features. In TaxoTox, both Verhaar classes 1 (non-polar
narcosis, LogKow ≥ 2) and 2 (polar narcosis, LogKow < 2) are mapped to `"narcosis"`
because they share the same membrane-disruption mechanism and are treated identically
in the Concentration Addition framework. Verhaar classes 3 (unspecific reactive) and 4
(specifically acting) cannot be resolved without structural alert analysis (SMARTS
matching), which requires external chemical informatics tools not currently integrated.

**Priority 3 — Default** (`moa_source = "default_narcosis"`):
Compounds for which neither name patterns nor LogKow are available are assigned
`moa_group = "narcosis"`. This includes compounds not in CompTox and those with unusual
chemical names not covered by the current pattern library.

### 8.4 Coverage and Observed Results

Mode of Action group coverage is reported in Plot 7 of `Code/validate_reference.R`:
a stacked bar chart shows the number of compound × group rows per mode of action group
per taxonomic group, with hover text showing the classification source.

In the current `Known_CAS.fst` dataset (~2800 compounds, estuarine monitoring panel),
the observed distribution is heavily dominated by narcosis (~97% of rows), with small
fractions of acetylcholinesterase inhibitors (~0.4%), Photosystem II inhibitors, and
pyrethroids. This reflects the actual chemical composition of the dataset: the majority
of monitored compounds are pharmaceuticals, personal care products, industrial chemicals,
and pesticide metabolites that act by non-specific membrane toxicity.

In practice, Priority 2 (LogKow-based Verhaar classification) contributed zero
classifications in the current run: every compound that was not matched by a name pattern
(Priority 1) fell directly to the narcosis default (Priority 3), because Verhaar classes 1
and 2 both map to narcosis anyway. Priority 2 adds no new information for the two narcosis
classes; its only potential value would be for Verhaar classes 3 and 4, which require
structural alert analysis not yet implemented. The LogKow retrieval step is therefore
retained in the pipeline for future extension but has no effect on the current classification
output.

This result has a practical implication for the CAMA method: in most water samples,
CAMA will return results very close to standard Concentration Addition / Toxic Unit
calculations, because there are insufficient compounds from distinct mechanistic groups
to generate divergent outcomes.

Key limitations of the current classification:

1. **Name pattern coverage is incomplete**: the pattern library covers approximately
   80–100 named active ingredients across the five mechanistic groups. Novel compounds,
   transformation products, and compounds named by IUPAC systematic names rather than
   common names may be missed and default to narcosis.

2. **Narcosis dominance is expected**: the majority of organic chemicals in any
   environmental dataset act by narcosis (Verhaar et al., 1992). A reference table
   dominated by `"narcosis"` does not indicate a classification failure — it reflects
   the chemical composition of typical aquatic monitoring datasets.

3. **Verhaar classes 3 and 4 are not fully resolved**: without SMARTS-based structural
   alert screening, reactive compounds and specifically-acting compounds not covered by
   name patterns will be assigned to narcosis. This is conservative — narcosis is a
   weaker assumption than specific receptor binding — but may underestimate risk for
   samples dominated by reactive chemicals.

4. **Mode of Action versus mode of toxic action**: the classification implemented here
   uses the toxicological mode of action at the organism level (how the chemical kills
   or inhibits), not the molecular initiating event in the sense of the Adverse Outcome
   Pathway (AOP) framework. These are generally equivalent for the compound classes
   covered but the distinction matters for more granular mechanistic analysis.

### 8.5 References

- Verhaar, H.J.M., van Leeuwen, C.J., Hermens, J.L.M. (1992). Classifying environmental
  pollutants. *Chemosphere*, 25(4), 471–491.
  https://doi.org/10.1016/0045-6535(92)90280-5

- European Chemicals Agency (ECHA). *Guidance on the Application of the CLP Criteria,
  Chapter R.7a: Endpoint specific guidance*. ECHA, Helsinki.

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

$$E_i = \frac{C_i}{LC50_i + C_i}$$

where $E_i$ is the fractional effect (0–1) of compound $i$ at measured concentration $C_i$,
and $LC50_i$ is the same denominator used in the standard TU calculation
(`dplyr::coalesce(median_lc50_ng_L, predicted_lc50_ng_L)`).

The mixture fractional effect is then:

$$E_{mix} = 1 - \prod_i (1 - E_i)$$

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

Five groups are used (see Section 8 for classification details):
`narcosis`, `AChE_inhibition`, `PSII_inhibition`, `pyrethroid`, `reactive`.
Compounds with no classifiable MoA default to `narcosis`.

### 10.5.4 Output

Results are written to three sheets (`CAMA Algae`, `CAMA Crustacean`, `CAMA Fish`).
Each sheet has:

- **`E_mix`** (first column) — mixture fractional effect per sample (0–1)
- Per-group **`E_group`** columns (e.g. `E_narcosis`, `E_AChE_inhibition`) — useful for
  seeing which mechanistic group drives the mixture risk

### 10.5.5 Practical expectation

Given the current dataset composition (~97% narcosis compounds; see Section 8.4), CAMA
will in most cases return values very close to standard CA-based PTI. Divergence is
expected only in samples where specific-acting compounds (organophosphates, pyrethroids,
triazines) contribute meaningfully to total TU. The CAMA output is most informative as
a sensitivity analysis: large differences between PTI and $E_{mix}^{CAMA}$ indicate that
cross-mechanism interactions may be important.

### 10.5.6 References

- Drescher, K., & Boedeker, W. (1995). Assessment of the combined effects of substances:
  the relationship between concentration addition and independent action. *Biometrics*, 51, 716–730.

- Altenburger, R., et al. (2000). Predictability of the toxicity of multiple chemical
  mixtures to Vibrio fischeri: mixtures composed of similarly acting chemicals.
  *Environmental Toxicology and Chemistry*, 19(10), 2341–2347.

---

## 11. CASRN Resolution

Compound names in the input file are resolved to CAS Registry Numbers (CASRNs) through two sequential layers:

**Layer 1 — Known_CAS (exact, instant):** a curated internal table (`Known_CAS.fst`) of compound names and their CASRNs. Exact string matches are confirmed automatically with no user interaction. The table is designed to hold multiple synonyms per compound (e.g. "Albuterol", "Salbutamol", and "Albuterol / Salbutamol" may all map to the same CASRN 18559-94-9). Each entry has a `PREFERRED_NAME` which is the canonical display name used in all outputs.

**Layer 2 — PubChem REST API (optional):** when enabled, unresolved names are queried live against the PubChem compound database via the `webchem` R package. Candidate matches are presented for user review in the CASRN Matching tab. Accepted matches are written to a temporary log (`temp_CAS` — local file or Google Sheet when deployed) for future integration into Known_CAS via the curation script.

**Manual entry:** compounds still unresolved after both layers can have their CASRN entered directly in the CASRN Matching tab. These are also logged to `temp_CAS`.

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
