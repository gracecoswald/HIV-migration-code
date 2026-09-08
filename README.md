# Human migration and the global distribution of HIV-1 genetic variants

Analysis code for *"Human migration and the global distribution of HIV-1 genetic variants"*

R scripts that link country-level HIV-1 variant distributions, derived from a
global molecular epidemiology database, to bilateral migration flow estimates,
and models the relationship between migration and variant-distribution
similarity.

---

## Quick start

Three input datasets cannot be shared here (see **Data** below). With the full
inputs in the working directory the whole analysis runs with:

```r
source("run_all.R")
```

With what is provided here, plus `migration_df_final.csv` (sent separately due to file size),
**10 of the 13 analysis scripts run**, reproducing Table 2, Appendix Table 4,
Supp Table 3 and all of the figures:

```r
source("00_libraries.R")
# copy `raw data/` and `produced datasets/` into the working directory
source("migration_RMSD_creation.R")
source("model_building.R")
source("country_diversity.R")
source("region_migration.R")
source("pie_chart_generation.R")
source("global_chord_diagrams.R")
source("africa_figures.R")
source("region_maps.R")
```

---

## 1. System requirements

### Software

| | Version |
|---|---|
| R | 4.5.3 |
| Operating system | Windows 11 x64 (build 26100) |

No other software is required. The code is plain R and should run on macOS and
Linux, but has only been tested on the platform above.

### R package dependencies

| Package | Version | | Package | Version |
|---|---|---|---|---|
| dplyr | 1.2.0 | | ggplot2 | 4.0.2 |
| tidyr | 1.3.2 | | RColorBrewer | 1.1-3 |
| tibble | 3.3.1 | | scales | 1.4.0 |
| stringr | 1.6.0 | | gridExtra | 2.3 |
| purrr | 1.2.1 | | patchwork | 1.3.2 |
| readxl | 1.4.5 | | sf | 1.0-24 |
| readr | 2.1.6 | | rnaturalearth | 1.2.0 |
| fixest | 0.14.2 | | rnaturalearthdata | 1.0.0 |
| alpaca | 0.3.5 | | png | 0.1-8 |
| modelsummary | 2.6.0 | | circlize | 0.4.18 |
| flextable | 0.9.10 | | vegan | 2.7-5 |
| corrplot | 0.95 | | grid | base (4.5.3) |

**`fixest` version matters.** Two coefficients in the full interaction models
sit in a badly conditioned design, and different `fixest` releases resolve it
differently. Version 0.14.2 reproduces the published standard errors; 0.13.2
does not. Pin this version if the exact values matter.

### Non-standard hardware

None. The analysis runs on a standard desktop or laptop; peak memory use is
under 4 GB.

### Other requirements

- `sf` requires the GDAL, GEOS and PROJ system libraries. These ship with the
  Windows and macOS CRAN binaries; on Linux they must be installed separately.
- `rnaturalearth` downloads map data on first use, so the figure scripts need
  an internet connection the first time they run.

---

## 2. Installation guide

No installation beyond R and the packages above.

```r
# from the repository directory
source("00_libraries.R")   # installs anything missing, then loads everything
```

**Typical install time on a normal desktop computer: 10-20 minutes** from a
clean R installation, downloading CRAN binaries. Most of that is `sf` and its
spatial dependencies. If the packages are already present, loading takes a few
seconds.

To pin the tested versions instead:

```r
install.packages("remotes")
remotes::install_version("fixest", "0.14.2")
```

---

## 3. Demo

The database export cannot be shared here, so the demo starts from the intermediate
datasets in `produced datasets/`, which were generated from the full dataset.

### Instructions to run

```r
# 1. Copy the contents of `raw data/` and `produced datasets/` into the
#    working directory, together with migration_df_final.csv (sent separately)
# 2. Then:
source("00_libraries.R")
source("migration_RMSD_creation.R")   # regenerates migration_RMSD.csv
source("model_building.R")            # Table 2, Appendix Table 4, Figure 4
source("country_diversity.R")         # diversity over time
source("region_migration.R")          # Supp Table 3
source("pie_chart_generation.R")      # pie charts and legends
source("global_chord_diagrams.R")     # Figure 1B
source("africa_figures.R")            # Figure 3
source("region_maps.R")               # Figure 2 and equivalents
```

### Expected output

| File | Contents |
|---|---|
| `migration_RMSD.csv` | Migration flows joined to RMSD; should match the copy in `produced datasets/` |
| `table2_main_results.docx` | Table 2 of the manuscript |
| `table4_standardised.docx` | Appendix Table 4 |
| `predicted_RMSD_prev.png` | Figure 4 |
| `complete_1995_diversity.csv` | Diversity indices over time with Kendall trend tests |
| `total_region_migration_new.csv`, `2015_mig_matrix.csv` | Supp Table 3 |
| `pie_charts_country/`, `pie_charts_region/`, legend PDFs | Variant pie charts |
| `chord_diagrams/` | Figure 1B and equivalents |
| `africa_maps/`, `chord_africa/` | Figure 3 and equivalents |
| `regional_maps/` | Figure 2 and equivalents |

Because `produced datasets/` was generated from the full dataset,
`model_building.R` **reproduces the published Table 2, Appendix Table 4 and
Figure 4 exactly**, given `fixest` 0.14.2.

`migration_RMSD_creation.R` regenerates `migration_RMSD.csv` from
`migration_df_final.csv` and `RMSDcountries.csv`; comparing it against the copy
in `produced datasets/` checks that step end to end.

### Expected run time

**Roughly 6-10 minutes** in total on a normal desktop computer, once packages
are installed. The pie chart generation accounts for most of it, every other
script runs in under two minutes.

### Scripts that cannot be run from what is provided

| Script | Needs |
|---|---|
| `subtyped_new.R` | the database export |
| `region_distributions.R` | `HIV_samples`, from `subtyped_new.R` |
| `migration_RMSD_adjusted.R` | `PLHIV_data.xlsx` |

Their outputs are supplied in `produced datasets/`, so the scripts downstream
of them still run.

---

## 4. Instructions for use

### Running on your own data

Supply your own versions of the four input files listed under **Data**,
keeping the same structure. The scripts read them by filename from the working directory, so
either match the filenames or edit the `read.csv` / `read_excel` line at the
top of the relevant script.

| File | Required structure |
|---|---|
| Database export | One row per study record, with `Site 1: Country` through `Site 30: Country`, `Year study started`, `Year study ended`, `Total number genotyped`, the `HIV-1 group M:` variant columns, the `CRF*` columns, the `URF* number` columns, and `Summary score for overall risk of study bias`. The genotyping fragment columns must sit at positions 343-352. |
| `migration_df_final.csv` | One row per origin-destination-year, with `orig`, `dest`, `orig_country`, `dest_country`, `year0`, `plot_area_origin`, `plot_area_dest`, and the flow estimates `da_min_open`, `da_min_closed`, `da_pb_closed`. |
| `population_data.csv` | World Bank format: `Country.Name`, `Series.Name`, `Series.Code`, and one column per year. |
| `PLHIV_data.xlsx` | UNAIDS format: country in column 4, one column per year. |

Country names must match across files, or pairs are silently dropped.
`migration_RMSD_creation.R` stops with an error naming any country that has no
region assigned; `migration_RMSD_adjusted.R` warns about countries with no
population match.

### Reproducing the manuscript

The database export and the PLHIV data cannot be shared, so the pipeline
cannot be run end to end from this repository alone. Everything downstream of
them can be, because their outputs are supplied in `produced datasets/`.

`model_building.R` run against those datasets reproduces **Table 2, Appendix
Table 4 and Figure 4 exactly**, given `fixest` 0.14.2. With
`migration_df_final.csv`, Supp Table 3 and all of the figures are regenerated
too, and `migration_RMSD.csv` can be regenerated and compared against the
provided copy as a check on that step.

Reproducing the variant distributions and the RMSD calculation themselves
requires the database export; the population-weighted flows require the PLHIV
data. Requests should be directed to the corresponding author.

---

## Data

### `/raw data/`

The scripts read these by filename from the working directory. Only
`population_data.csv` is included in this repository.

| File | Contents | Provided? |
|---|---|---|
| `HIVIDDONDPHProject_DATA_LABELS_2024-05-07_1226.xlsx` | Export of the global molecular epidemiology database, one row per study record. | **No.** Requests to the corresponding author. |
| `migration_df_final.csv` | Bilateral migration flow estimates by country pair and time period. | **Sent separately**, too large for GitHub. |
| `population_data.csv` | World Bank country populations by year. | **Yes**, in `raw data/`. |
| `PLHIV_data.xlsx` | UNAIDS people living with HIV, by country and year. | **No.** Not ours to redistribute. Requests to the corresponding author. |

### `/produced datasets/` — provided

The intermediate datasets produced by running the pipeline on the full raw
dataset. These **are** included, so the output of each stage can be inspected,
and so the modelling step can be run and checked without the raw inputs.

| File | Written by |
|---|---|
| `counts_country_subtype_distirbutions.csv` | `subtyped_new.R` |
| `proportions.csv` | `subtyped_new.R` |
| `RMSDcountries.csv` | `subtyped_new.R` |
| `indices.csv` | `subtyped_new.R` |
| `counts_region_time.csv` | `region_distributions.R` |
| `proportions_region_time.csv` | `region_distributions.R` |
| `migration_RMSD.csv` | `migration_RMSD_creation.R` |
| `migration_adj_RMSD.csv` | `migration_RMSD_adjusted.R` |
| `migration_adj_RMSD_absdiff.csv` | `migration_RMSD_adjusted.R` |

Copy these into the working directory to run the scripts listed under
**Demo** above.

`net_flow_adj.csv`, also written by `migration_RMSD_adjusted.R`, is not
included, nothing downstream reads it.

---

## Scripts

Run in this order. `run_all.R` does it for you.

### Setup

| Script | Purpose |
|---|---|
| `00_libraries.R` | Loads every package used anywhere. No analysis script contains `library()` calls. Sourced by `run_all.R`. |
| `ne_country_names.R` | Maps our country names to `rnaturalearth`'s (`"Congo (the Democratic Republic of the)"` to `"Dem. Rep. Congo"`). Sourced by `region_maps.R`, not called directly. |

### Core analysis

| Script | Reads | Writes | Purpose | Run time |
|---|---|---|---|---|
| `subtyped_new.R` | database xlsx | `counts_country_subtype_distirbutions.csv`, `proportions.csv`, `RMSDcountries.csv`, `indices.csv` | Cleans the database, assigns regions, collapses variant columns, splits multi-year studies across years, aggregates by country and time period, computes pairwise RMSD between every country pair, and computes Shannon and Simpson diversity indices. | 10-15 min |
| `migration_RMSD_creation.R` | `migration_df_final.csv`, `RMSDcountries.csv` | `migration_RMSD.csv` | Sums bidirectional migration per country pair, assigns regions to both sides, joins to RMSD. Also builds `net_flow` and `final_region_flow`, used by the figure scripts. | < 1 min |
| `migration_RMSD_adjusted.R` | `migration_df_final.csv`, `population_data.csv`, `PLHIV_data.xlsx` | `net_flow_adj.csv`, `migration_adj_RMSD.csv`, `migration_adj_RMSD_absdiff.csv` | Weights flows by destination population (per 100,000), and computes the absolute difference in HIV prevalence between each pair. | < 1 min |
| `model_building.R` | `migration_RMSD.csv`, `migration_adj_RMSD.csv`, `migration_adj_RMSD_absdiff.csv` | `table2_main_results.docx`, `table4_standardised.docx`, `predicted_RMSD_prev.png` | Fixed-effects quasi-Poisson models of RMSD on migration flow, prevalence difference and time. Produces Table 2, Appendix Table 4 and Figure 4. | < 1 min |

### Descriptive tables

| Script | Writes | Purpose | Run time |
|---|---|---|---|
| `country_diversity.R` | `div_mig_popn.csv`, `complete_1995_diversity.csv` | Diversity indices over time for countries with data in every period, with Kendall trend tests. | < 1 min |
| `region_distributions.R` | `counts_region_time.csv`, `proportions_region_time.csv` | Regional variant distributions by time period (Supp Table 2). | < 1 min |
| `region_migration.R` | `total_region_migration_new.csv`, `2015_mig_matrix.csv` | Immigrant / emigrant / within-region totals per region and period, and the region-to-region flow matrix for 2015-2019 (Supp Table 3). Matrices for the other periods are held in the `mats` list. | < 1 min |

### Figures

| Script | Writes | Purpose | Run time |
|---|---|---|---|
| `pie_chart_generation.R` | `pie_charts_country/`, `pie_charts_region/`, `Subtype_Legend.pdf`, `Small_Subtype_Legend.pdf` | One variant pie chart per country and per region per period, plus the colour legends. Builds `country_pie_files`, used by the map scripts. | 3-5 min |
| `global_chord_diagrams.R` | `chord_diagrams/` | Global chord diagrams of regional migration flows, one per period plus a pooled version (Figure 1B). | < 1 min |
| `africa_figures.R` | `africa_maps/`, `chord_africa/` | Africa maps with migration flows and country pie charts, and Africa chord diagrams, one of each per period (Figure 3). | < 1 min |
| `region_maps.R` | `regional_maps/` | Maps for the Americas, Europe and Africa, Europe and North America, Europe, and Asia, one per period (Figure 2 and equivalents). | 1-2 min |

The figure scripts must run in this order: `pie_chart_generation.R` builds the
pie charts that `africa_figures.R` and `region_maps.R` overlay on the maps.

### Sensitivity analyses

Each is a self-contained rerun of the entire pipeline with one filter applied,
writing outputs suffixed `_pol` or `_rob`. They mirror the main scripts exactly
apart from that filter. **Run in a clean R session** — they redefine objects
used by the main analysis.

| Script | Filter |
|---|---|
| `sensitivity_pol.R` | Only records where every genotyping fragment is *pol* or full length |
| `sensitivity_rob.R` | Excludes records assessed as high risk of bias |

---

## Output naming

Main analysis outputs are unsuffixed. Sensitivity outputs carry `_pol` and
`_rob`, so they never overwrite the main results.
