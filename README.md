# backcalc

> **Archived.** No longer maintained; last updated June 2020. It depends on the retired `EpiNow` and `NCoVUtils` packages.

Evaluation of approaches to back-calculating case counts by date of infection from case counts by date of report. The report compares a convolution-based approach in Stan with the independent sampling method then used in `EpiNow`, on simulated data and on reported COVID-19 cases from four countries, and adds an early in-model estimate of Rt using a Gaussian process.

- `report.Rmd`: the analysis; `report.md` and `report.html` are its rendered versions, with figures in `figures/`.
- `nowcast.stan`, `nowcast.R`: the Stan model and the R code that fits it.
- `R/generate_simulations.R`, `scripts/generate-outbreak-scenarios.R`: generate the simulated outbreaks.
- `data/`: simulated scenarios and case counts.
