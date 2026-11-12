# HRP2 persistence: publication code

Code for *Estimating HRP2 antigen persistence and effect on rapid diagnostic test data: a Bayesian hierarchical approach* (Walters et al.).

The analysis uses paired rapid diagnostic test (RDT) and light microscopy (LM) results from Demographic and Health Surveys (DHS) and Malaria Indicator Surveys (MIS) in sub-Saharan Africa. It has two parts:

- **Regression analysis**: a Bayesian categorical mixed-effects regression of the joint RDT and LM result, giving the false positive rate (FPR) and false discovery rate (FDR) by covariate.
- **Calibration analysis**: *malariasimulation* run over a grid of entomological inoculation rate (EIR) and HRP2 decay, Gaussian process (GP) emulators of each admin-1 region's likelihood surface, and a Bayesian hierarchical model of HRP2 half-life across regions.

The repository is self-contained: every figure, table and number in the manuscript is written by code in `src/`, except the hand-drawn *malariasimulation* model framework (`figures/static/`). Survey microdata and site files are restricted and must be obtained separately (see [Data access](#data-access)).

## Repository structure

```
mal_proj_publication/
├── launch.R                  runs any pipeline stage (Rscript launch.R list)
├── run_sections.R            runs one step script with settings overrides
├── provision.R, pkgdepends.txt   cluster package provisioning (hipercow)
├── .Renviron.example         template for your DHS login and optional paths
├── src/
│   ├── 1_data_prep/          DATA CLEANING AND TRANSFORMATION
│   ├── 2_regression/         REGRESSION ANALYSIS
│   ├── 3a_calibration_grid/  CALIBRATION ANALYSIS: simulation grid and likelihood
│   ├── 3b_calibration_emulator/  CALIBRATION ANALYSIS: GP emulators
│   ├── 3c_calibration_mcmc/  CALIBRATION ANALYSIS: hierarchical MCMC
│   │   └── sensitivity/      SENSITIVITY ANALYSES (also in 3a/ and 3b/)
│   ├── 4_sample_posterior/   CALIBRATION ANALYSIS: posterior predictive check
│   ├── plots.R               manuscript figures -> figures/
│   ├── results_numbers.Rmd   every cited number and generated table -> results/numbers.tex
│   ├── export_posteriors.R   posterior draws and summaries -> results/posteriors/
│   ├── sync_manuscript.R     copies figures and numbers.tex into your manuscript folder
│   └── functions/            function libraries
├── data/                     public inputs (tracked); data/site_files/ is restricted (ignored)
├── figures/                  code-generated figures; figures/static/ holds the hand-drawn model framework
├── results/                  numbers.tex, LOO table, posteriors/
```

## Analysis components

Every script's header states its component on an `ANALYSIS COMPONENT:` line.

| Component | Scripts |
|---|---|
| **Data cleaning and transformation** | `src/1_data_prep/01_format_dhs_mis_data.R` (DHS/MIS download and indicators), `02_survey_design_effects.R`, `03_structure_mis_data.R` (site matching, deletion risk, seasonality); `src/functions/central_functions.R`, `check_fns.R` |
| **Regression analysis** | `src/2_regression/fit_fdr_fpr_regression.R` (main model, FPR/FDR, average marginal effects), `src/2_regression/fixed_effect_selection.R` (LOO); `src/functions/brms_int_cluster_fns.R`, `fpr_fdr_map_fns.R` |
| **Calibration analysis** | `src/3a_calibration_grid/grid_simulation.R`, `src/3b_calibration_emulator/emulation.R`, `src/3c_calibration_mcmc/hierarchical_mcmc.R`, `src/4_sample_posterior/sample_posterior.R`; `src/functions/cluster_fns.R`, `emulator_fns.R`, `validate_emulators.R`, `plot_emulators.R`, `hierarchical_mcmc_fns.R`, `posterior_summary_fns.R`, `plot_mcmc_diagnostics.R` |
| **Sensitivity analyses** (calibration) | RDT sensitivity 0.8 and 1.0 (main 0.9): `3a_calibration_grid/sensitivity/rdt_sensitivity_log_lik.R`, `3b_calibration_emulator/sensitivity/rdt_sensitivity_emulation.R`, `3c_calibration_mcmc/sensitivity/rdt_sensitivity_mcmc.R`. Prior SD of sigma_P 0.05 and 0.1 (main 0.01): `3c_calibration_mcmc/sensitivity/sigma_prior_sd_mcmc.R`. These re-run the step script with different settings and hold no model code of their own. |
| **Outputs** | `src/plots.R`, `src/results_numbers.Rmd`, `src/export_posteriors.R`, `src/sync_manuscript.R`; `src/functions/numbers_fns.R`, `tex_table_fns.R`, `schematic_fns.R`, `data_stream_fns.R` |
| **Pipeline control** | `launch.R`, `run_sections.R`, `src/functions/pipeline_fns.R`, `paths.R`, `provision.R`, `pkgdepends.txt` |

## How the code maps to the manuscript

### Main text

| Manuscript item | Output file | Produced by |
|---|---|---|
| Figure 1, HRP2 dynamics schematic (`fig:dynamics`) | `figures/schematic.png` | `plots.R`, `PLOT_SCHEMATIC` (`schematic_fns.R`; all three panels share one simulated cohort) |
| Figure 2, FPR and FDR maps (`fig:regression`) | `figures/fpr_fdr_maps_ssa.png` | `plots.R`, `PLOT_FPR_FDR_MAPS` (from the step 2 fit) |
| Figure 3, RDT:LM ratio and half-life forest (`fig:halflife`) | `figures/combined_ratio_density_and_forest_plot.png` | `plots.R`, `PLOT_DENSITY_FOREST` (from the step 3c main fit) |
| Figure 4, guidance for programmes (`fig:NMP_guidance`) | manuscript `figures/nmcp_guidance_body.tex` (TikZ, kept with the manuscript) | numbers from `results_numbers.Rmd` (`res.*`, `fig.*` keys) |
| Every number in the Methods and Results | `results/numbers.tex` | `results_numbers.Rmd`; cited as `\nm{key}` |

### Appendix

| Appendix section | Item | Produced by |
|---|---|---|
| Data: household survey data | Data stream figure (`fig:data_stream`) and counts | `plots.R`, `PLOT_DATA_STREAM`; counts from `data_stream_fns.R`, shared with `results_numbers.Rmd` (`data.*` keys) |
| Data: *pfhrp2* gene deletion | Deletion table (`tab:deletion`) | inputs `data/hrp2_deletion.csv`, `data/pfhrp2_deletion_risk.csv`; table currently typed in `appendix.tex` |
| Regression: fixed effect selection | ELPD table (`tab:elpd`) and text | `fixed_effect_selection.R` -> `results/table_loo_fixed_effect_selection.csv` -> `results_numbers.Rmd` (`tab.elpd`, `elpd.*`) |
| Regression: average marginal effects | `fig:app_fpr_fdr` | `plots.R`, `PLOT_AME` -> `figures/fpfn_ame_forest.png` |
| Calibration: *malariasimulation* | Model framework (`fig:malsim_schematic`) | hand-drawn: `figures/static/model_framework.pptx` -> `.png` |
| Calibration: GP emulators | Emulator flags table (`tab.emuflags`) | `emulation.R` -> `results_numbers.Rmd` |
| Calibration: hierarchical model | MCMC evaluation text and table (`tab.mcmc`) | `hierarchical_mcmc.R` -> `results_numbers.Rmd` (`mcmc.*`) |
| Sensitivity analysis | sigma_P prior table (`tab.sigma`) and text | `sigma_prior_sd_mcmc.R` -> `results_numbers.Rmd` (`post.*.sd005`, `post.*.sd010`) |
| Posterior predictive coverage | `fig:scatter` | `sample_posterior.R` -> `plots.R`, `PLOT_PREVALENCE_WU_SCATTER` -> `figures/prevalence_rdt_lm_wu_scatter.png` |
| Population and site results | `tab.post`, `tab.sites` | `hierarchical_mcmc.R` -> `results_numbers.Rmd` |

When `HRP2_MANUSCRIPT_DIR` points to a local copy of the manuscript, `results/numbers.tex` holds only the keys its `.tex` files cite, and `results_numbers.Rmd` reports any cited key the code does not produce. Without it, `numbers.tex` holds every computed key.

### Posterior estimates

`src/export_posteriors.R` writes every posterior estimate to `results/posteriors/`:

- `regression/`: every draw and a summary (mean, SD, median, 95% credible interval, R-hat, effective sample size) of each fixed effect, random-effect SD and survey random effect, plus the FPR, FDR and average marginal effects.
- `calibration/<fit>/`: for the main fit and each sensitivity fit, every pooled draw of mu_P, sigma_P and the population median half-life; thinned joint draws of each region's half-life and EIR; region, country and unobserved-region summaries; convergence diagnostics; priors.

`<fit>` is `main_rdt_se_090`, `sensitivity_rdt_se_080`, `sensitivity_rdt_se_100`, `sensitivity_sigma_sd_0.05` or `sensitivity_sigma_sd_0.1`. These files hold model parameters and model-derived quantities only, never survey microdata or site-file inputs.

## Data access

Two inputs are restricted and are **not** in this repository. `.gitignore` excludes them and everything derived from them (`outputs/`, `data/site_files/`, `rdhs.json`, the raw DHS cache and raw DHS file formats).

### DHS/MIS microdata

1. Register with the DHS Program: https://dhsprogram.com/data/new-user-registration.cfm
2. Request access for your own project to the surveys in sub-Saharan Africa with malaria testing (household member, children's and biomarker recode files).
3. Copy `.Renviron.example` to `.Renviron` in the project root and add your login:

   ```
   DHS_EMAIL=you@institution.org
   DHS_PROJECT=Your DHS project title, exactly as registered
   ```

4. Restart R. Step 1a (`01_format_dhs_mis_data.R` with `REPULL_DATA = TRUE`) asks for your DHS password and downloads the data with `rdhs`. Your password is never written to the repository.

### Site files

The *malariasimulation* site files (one `<ISO>_site.rds` per country, and the Uganda subnational site file `UGA_snt_site.RDS`) are available on request from the corresponding author. Most country site files can also be built with the `site` package (`remotes::install_github("mrc-ide/site@site-2601")`, then `site::fetch_site()`), which needs access to the MRC Centre for Global Infectious Disease Analysis site-file store.

Put the files in `data/site_files/`, or point to them with `HRP2_SITE_CACHE_DIR` and `HRP2_UGA_SNT_SITE_FILE` in `.Renviron`.

### Public inputs (tracked)

| File | Source |
|---|---|
| `data/pfhrp2_deletion_risk.csv` | Modelled *pfhrp2* deletion risk, Watson et al. 2025 |
| `data/hrp2_deletion.csv` | Pooled *pfhrp2* deletion prevalence, Thomson et al. 2020 systematic review |
| `data/name_1_map.csv` | Uganda admin-1 name matching |
| `data/geoboundaries/` | Admin boundaries from geoBoundaries (CC BY 4.0) |

## Setup

1. Install R and the packages below.
2. Copy `.Renviron.example` to `.Renviron` and fill it in (see [Data access](#data-access)).
3. Open `hrp2_persistence.Rproj` (or set the working directory to the project root). Every script assumes the working directory is the project root.

**Packages.** CRAN: `data.table`, `ggplot2`, `patchwork`, `ggtext`, `ggrepel`, `ltc`, `scales`, `sf`, `rnaturalearth`, `rnaturalearthhires`, `jsonlite`, `dplyr`, `stringr`, `stringdist`, `tidyverse`, `haven`, `survey`, `rdhs`, `brms`, `rstan`, `loo`, `posterior`, `DiceKriging`, `individual`, `knitr`, `rmarkdown`, `broom.mixed`. GitHub: `mwalte10/malariasimulation@age_spec_states` (the *malariasimulation* fork used for the grid), `mrc-ide/netz`, `mrc-ide/site@site-2601`, `mrc-ide/hipercow`. `cmdstanr` is optional (the cluster brms backend; see `pkgdepends.txt`).

## Running the analysis

`launch.R` runs every stage. From the project root:

```
Rscript launch.R list                 # every action and what it does
Rscript launch.R 1_data_prep          # run one or more actions in order
```

Or interactively: `RUN <- "plots"; source("launch.R")`.

| Order | Actions | Where |
|---|---|---|
| 1 | `1_data_prep` | local |
| 2 | `2_regression_fit`, `2_regression_selection` | local (hours per model) |
| 3 | `3a_cache_submit`, `3a_cache_collect`, `3a_grid_submit`, `3a_grid_collect`, `3a_likelihood` | DIDE cluster |
| 4 | `3b_submit`, `3b_collect` | DIDE cluster |
| 5 | `3c_run`, `3c_collect`, `3c_pool` | DIDE cluster |
| 6 | `4_submit`, `4_collect`, `4_aggregate` | DIDE cluster |
| 7 | `plots`, `results_numbers`, `export_posteriors`, `manuscript` | local (minutes) |

Step 7 rebuilds everything the manuscript uses from the saved fits, so it can be re-run on its own once steps 1 to 6 have finished.

**Cluster runs.** Steps 3a to 4 submit tasks with `hipercow`. Cluster tasks read the `.Renviron` in the project root, so keep `outputs/` and `data/site_files/` inside the project (the defaults) or use network paths that the cluster nodes can reach.

## What git tracks

| Tracked | Ignored |
|---|---|
| `src/`, `launch.R`, `run_sections.R`, provisioning files | `outputs/` (all DHS/MIS-derived data, simulations, emulators, fits, MCMC chains) |
| `data/` public inputs | `data/site_files/` and any `*_site.rds` |
| `figures/` (code-generated) and `figures/static/` | `.Renviron`, `rdhs.json`, raw DHS files and cache |
| `results/` (numbers.tex, LOO table, posteriors) | `archive/`, `logs/`, `hipercow/` |
| | `manuscript/` (kept outside the repository), knitted HTML |

## Building the manuscript

The LaTeX sources (`main.tex`, `appendix.tex`, the guidance TikZ figure and the bibliography) are kept outside this repository. To update them from the code, set `HRP2_MANUSCRIPT_DIR` in `.Renviron` to your local copy (for example an Overleaf git clone) and run `Rscript launch.R manuscript`, which copies `results/numbers.tex` and the figures into it. The manuscript cites each number as `\nm{key}`; a key that `numbers.tex` does not define prints as a red `??key`.
