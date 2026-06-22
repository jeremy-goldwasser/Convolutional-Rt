# code/

Analysis code for *Fast, Frequentist Estimation of Epidemic Reproduction Numbers*.

## Layout

```
code/
├── _paths.R               Auto-detects the repo root from any working directory;
│                          defines data_dir, results_dir, figures_dir.
├── install_packages.R     Installs all CRAN dependencies.
├── examples/
│   └── demo.R             Self-contained worked example: simulate an outbreak,
│                          fit ConvRt retrospectively and in real time, plot R(t).
└── analyses/              One subdirectory per dataset / comparison.
    ├── covidestim/        COVID-19 comparison against covidestim (CA, LA, counties).
    ├── flu/               U.S. flu hospitalization analysis (seasons 2022/23, 2023/24).
    ├── rtestim/           Benchmark against rtestim on renewal- and SEIR-simulated data.
    ├── sim_flu/           Flu-based SEIR simulation benchmark (retro + real-time).
    ├── tf/                Trend-filtering variant comparison.
    └── weekly/            Weekly-aggregated estimation proof-of-concept.
```

Each subdirectory contains a `_common.R` that loads ConvRt, detects the repo root, and defines dataset-specific helpers shared across that subdirectory's scripts. Analysis scripts are standalone: run them with `Rscript` from any directory and they write results to `results/` and figures to `figures/`.
