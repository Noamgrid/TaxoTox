# TaxoTox User Guide

TaxoTox is an R Shiny application for aquatic mixture toxicity screening. Given a table of
pollutant concentrations, it calculates Toxic Units (TU) and the Pollution Toxicity Index
(PTI) for algae, crustaceans, and fish, using species-sensitivity data from the US EPA
ECOTOX Knowledgebase.

This guide walks through using the app end to end.

**Live app**: [https://yairsuari.shinyapps.io/TaxoTox/](https://yairsuari.shinyapps.io/TaxoTox/)


---

## Quick Start

Three steps, in the sidebar, in order:

1. **1. Upload Samples Data** — select the correct data layout, then upload your file.
2. **2. Calculate Toxicity** — becomes active once all compounds have a CASRN.
3. **3. Download Results** — exports a multi-sheet Excel workbook.

The sections below cover each step in detail.

---

## Running TaxoTox

Open the hosted app — nothing to install: [https://yairsuari.shinyapps.io/TaxoTox/](https://yairsuari.shinyapps.io/TaxoTox/)

---

## Preparing Your Input File

**Accepted formats**: CSV, TSV, TXT, XLS, XLSX.

> ⚠️ **All concentrations must be in ng/L** (nanograms per litre). The app does not
> convert units for you — values entered in µg/L, mg/L, ppb, or any other unit will be
> silently treated as ng/L and produce incorrect Toxic Unit / PTI results. Convert your
> data to ng/L before uploading.

The app accepts two layouts and auto-detects which one you're using:

**Layout A — pollutant names in the first column (default)**

| Compound    | Station1 | Station2 | ... |
|---|---|---|---|
| Caffeine    | 10.5 | 18.2 | ... |
| Atrazine    | 2.1  | 0.5  | ... |
| Bisphenol A | 2.1  | 0.3  | ... |

**Layout B — station names in the first column**

| Station  | Caffeine | Atrazine | Bisphenol A | ... |
|---|---|---|---|---|
| Station1 | 10.5 | 2.1 | 2.1 | ... |
| Station2 | 18.2 | 0.5 | 0.3 | ... |

Layout B is automatically transposed to Layout A before processing. If more than 75% of
the names in the first column aren't recognised as compounds, the app will offer to
transpose the data for you — a useful signal that the wrong layout was selected.

---

## Step 1 — Upload Your Data

In the sidebar, choose the matching **Data layout** radio button, then click
**"1. Upload Samples Data"**. The app immediately resolves compound names to CAS Registry
Numbers (CASRNs) and navigates to the **CASRN Matching** tab.

> 🎥 **Screencast — Preparing & Uploading Your Data**
> Shows: formatting an input file (choosing Layout A vs Layout B, confirming concentrations
> are in ng/L), selecting the matching data layout in the sidebar, uploading the file, and
> the automatic jump to the CASRN Matching tab.
> `[ADD YOUTUBE LINK]`

---

## Step 2 — Resolve CASRN Matches

Compound names are resolved through up to three layers:

1. **Known_CAS (instant, always on)** — a curated internal table of commonly monitored
   compounds. Exact matches are confirmed automatically.
2. **PubChem REST API (optional)** — toggle with the **"Web search for pollutant names"**
   checkbox. Candidate matches appear in the **Candidate Matches** table; review them,
   uncheck any incorrect suggestions, then click **"Submit Selected Lines"**.
3. **Manual CASRN Entry** — for anything still unresolved, enter the CASRN directly
   (e.g. `107-13-1`) in the section below the table and click **"Add CASRN"**.

The **"2. Calculate Toxicity"** button activates once every compound has a confirmed CASRN.

If two rows resolve to the same CASRN (e.g. a compound listed under two synonyms), the app
shows a duplicate-CASRN warning — both rows will still be summed into the PTI, effectively
double-counting that compound, so remove one before proceeding.

> 🎥 **Screencast — Resolving CASRN Matches**
> Shows: reviewing PubChem candidates, submitting selected matches, and manually entering
> a CASRN for an unmatched compound.

[![Watch: Resolving CASRN Matches](https://img.youtube.com/vi/Md_Fhz-2pRQ/0.jpg)](https://youtu.be/Md_Fhz-2pRQ)

---

## Step 3 — Calculate Toxicity

Click **"2. Calculate Toxicity"**. For each compound *i*, sample *j*, and taxonomic group
*g* (algae, crustacean, fish):

```
TUᵢⱼ = Cᵢⱼ (ng/L) / median LC50ᵢ (ng/L)
PTIⱼ = Σᵢ TUᵢⱼ
```

**Risk interpretation:**

| PTI | Interpretation |
|---|---|
| ≥ 1.0 | Acute risk — mixture may cause lethal or immobilisation effects |
| 0.1 – 1.0 | Chronic risk — sub-lethal effects on sensitive species possible |
| < 0.1 | Low risk — mixture unlikely to cause measurable toxicity |

Compounds with no ECOTOX data but a modelled prediction (CompTox/OPERA) are flagged with
a `~` prefix (e.g. `~ibuprofen`) in table headers and Excel exports, denoting extra
uncertainty in that value.

**Where to find the results:**

- **In the app** — results populate immediately in the **Toxicity Plots** tab (top 10
  riskiest samples and top 10 most toxic pollutants per taxonomic group, interactive) and
  the **Toxicity Tables** tab (full per-compound TU/PTI tables, colour-coded by the risk
  thresholds above). See *Step 4 — Viewing Results* below for details.
- **In Excel** — clicking **"3. Download Results"** exports the same per-compound TU/PTI
  values (plus the original data, CASRN report, and toxicity coverage) as a multi-sheet
  workbook. See *Step 4 — Viewing Results → Download Results* below for the full sheet
  layout.

> 🎥 **Screencast — Calculating Toxicity**
> Shows: clicking "2. Calculate Toxicity", the results populating the Toxicity Plots and
> Toxicity Tables tabs, and the corresponding output in the downloaded Excel workbook.

[![Watch: Calculating Toxicity](https://img.youtube.com/vi/oRKlXdP2H1E/0.jpg)](https://youtu.be/oRKlXdP2H1E)

---

## Step 4 — Viewing Results

- **Toxicity Plots tab** — top 10 riskiest samples by PTI, and top 10 most toxic pollutants,
  per taxonomic group. All plots are interactive (hover for details, click legend entries
  to isolate groups).
- **Toxicity Tables tab** — full per-compound TU / PTI tables per group, colour-coded by
  the risk thresholds above (orange ≥ 0.1, red ≥ 1.0).

Read the disclaimer at the top of the Toxicity Plots tab before interpreting results —
TaxoTox is a research screening tool, not a basis for regulatory decisions on its own.

> 🎥 **Screencast — Viewing Results**
> Shows: navigating the Toxicity Plots and Toxicity Tables tabs, hovering a plot,
> filtering a table, then clicking "3. Download Results" and opening the resulting
> workbook (or zip, if Advanced Assessment was used).

[![Watch: Viewing Results](https://img.youtube.com/vi/KDTR450KU7w/0.jpg)](https://youtu.be/KDTR450KU7w)

### Download Results

Click **"3. Download Results"**. The standard export is a single Excel workbook with six
sheets:

- **Original data** — the uploaded concentration matrix (ng/L)
- **Algae / Crustacean / Fish toxicity** — per-compound TU and PTI per sample
- **CASRN Report** — matched and unmatched compounds with CASRNs
- **Toxicity Coverage** — which taxonomic groups have ECOTOX data for each compound

If any Advanced Assessment method is enabled, a second workbook (one sheet per method ×
taxonomic group) is produced and both files are bundled into a single `.zip`.

---

## Optional — Advanced Assessment

Expand **"Advanced Assessment"** in the sidebar to opt into additional methods, run
alongside the standard calculation on the same button press:

- **HC5-based TU/PTI** — uses the 5th-percentile Species Sensitivity Distribution
  (HC5) instead of the median LC50 as the denominator, a more protective, regulatory-style
  threshold. Produces systematically higher PTI values than the standard method.
- **Benchmark Hazard Index** — replaces the LC50 denominator with a national regulatory
  benchmark: **US EPA Aquatic Life Benchmarks**, the only framework this method offers —
  checking the box is all that's needed, there's no separate framework to pick. Three
  other frameworks (EU EQS, AU ANZG, CA CCME) were evaluated and removed because their
  compound coverage in TaxoTox's reference table was too low to be useful (each under 3%,
  vs. US EPA's ~9–10%) — see `TaxoTox_Technical_Methods.md` Section 6.3.
- **Independent Action (IA)** — an alternative to the default Concentration Addition
  model; treats compounds as acting independently rather than additively. More appropriate
  when a mixture is dominated by one or two potent, specifically-acting compounds.
- **CAMA** — applies Concentration Addition within groups of compounds sharing the same
  mode of action, then Independent Action across groups. Mode of Action groups come from
  a real, externally curated database (Kramer et al. 2024) rather than a fixed category
  list, and compounds it doesn't classify get an explicit `unknown` group rather than
  being assumed to act by baseline narcosis. Real classifications currently cover only
  ~11% of TaxoTox's compound universe, so CAMA produces **two output variants** per
  taxonomic group: one including the `unknown`-MoA compounds as their own group, and one
  ("known MoA") excluding them entirely — compare the two to see how much the
  unclassified compounds are driving the result. Expect CAMA's `E_mix` to track the
  standard PTI closely for most samples regardless of coverage — the two mathematically
  converge whenever every mode-of-action group's toxic-unit sum stays well below 1 (true
  for the great majority of real water samples); CAMA meaningfully diverges from PTI only
  when one mode-of-action group is concentrated enough on its own to approach or exceed a
  toxic-unit sum of 1. See `TaxoTox_Technical_Methods.md` Section 10.5 for details and a
  worked diagnostic against TaxoTox's own sample datasets.
- **Taxon-Sensitive PTI (Nowell et al. 2014)** — uses the USGS NAWQA Pesticide Toxicity
  Index (PTI) methodology (Nowell, Norman, Moran, Martin & Stone, 2014, *Sci. Total
  Environ.* 476–477:144–157), the reference method behind the published NWQN/RSQA PTI
  values many USGS water-quality studies report. The denominator ("sensitive toxicity
  concentration", STC) is read directly from Nowell et al.'s own published reference
  table (their Supplementary Appendix B), not recomputed from ECOTOX — so results
  exactly reproduce, rather than approximate, the values reported in downstream USGS
  literature (e.g. Covert et al. 2020). Covers **fish and cladocerans only** (Nowell's
  method doesn't include algae, and TaxoTox doesn't model their third taxon, benthic
  invertebrates); the cladoceran denominator is restricted to the 17 cladoceran genera
  Nowell et al. tested, not TaxoTox's full crustacean group, so results are directly
  comparable to published cladoceran-PTI values in a way the Benchmark Hazard Index
  method (a cross-invertebrate regulatory value) is not. Coverage is narrower than the
  other methods since it's limited to the ~480 fish and ~460 cladoceran-genus
  pesticides/degradates Nowell et al. published STC values for. See
  `TaxoTox_Technical_Methods.md` Section 10.6 for the full algorithm.

Full formulas, data sources, and validation for each method are in
`TaxoTox_Technical_Methods.md` (Sections 5, 6, 10, 10.5).

> 🎥 **Screencast — Advanced Assessment**
> Shows: enabling an advanced method (e.g. Benchmark Hazard Index) and locating the
> corresponding sheet in the downloaded workbook.

[![Watch: Advanced Assessment](https://img.youtube.com/vi/9EwMXNivGl0/0.jpg)](https://youtu.be/9EwMXNivGl0)

---

## Troubleshooting / FAQ

**Nothing matched after upload / most compounds are unresolved.**
Check that you selected the correct data layout (Layout A vs B). If more than 75% of the
first-column names are unrecognised, the app will prompt you to transpose automatically.

**The "2. Calculate Toxicity" button stays greyed out.**
At least one compound still needs a CASRN. Check the CASRN Matching tab for remaining
unmatched entries.

**I see a duplicate-CASRN warning.**
Two rows in your file resolved to the same compound (usually two synonyms). Remove one
before calculating, or the compound's concentration will be double-counted in the PTI.

**A compound in my results has a `~` in front of its name.**
No ECOTOX experimental data exists for that compound; the denominator is a modelled
(CompTox/OPERA) prediction instead, which carries additional uncertainty.

**A compound has no PTI contribution at all in one taxonomic group.**
Neither ECOTOX data nor a model prediction exists for that compound in that group (most
commonly algae, since no OPERA endpoint exists for algae). It contributes zero to that
group's PTI, which may underestimate risk — see the Limitations section of
`TaxoTox_Technical_Methods.md`.

**Can I contribute a new CASRN back to the shared reference table?**
Yes — manually entered CASRNs are logged automatically for later curation into
`Known_CAS`. No action needed on your part beyond entering them correctly.

---

## Citing TaxoTox

See the accompanying journal paper for the citation. Methodological detail, data
provenance, and validation results are documented in
[`TaxoTox_Technical_Methods.md`](TaxoTox_Technical_Methods.md).
