# DeadWood_Example

Example repository showcasing how to use these modules. Built
with [SpaDES](https://spades.predictiveecology.org/).

## Modules

Three modules are downloaded from the
[BosunForestEcology](https://github.com/BosunForestEcology) GitHub
organisation when `global.R` is run:

| Module | Description |
|--------|-------------|
| [DeadWood_snagDecay](https://github.com/BosunForestEcology/DeadWood_snagDecay) | Advances snags through DC1–DC5 via a Markov transition matrix; transfers fallen snags to the DWD pool |
| [DeadWood_DWDDecay](https://github.com/BosunForestEcology/DeadWood_DWDDecay) | Manages downed woody debris; receives fallen snags, advances DC, removes fully decomposed material |
| [DeadWood_Biomass](https://github.com/BosunForestEcology/DeadWood_Biomass) | Converts decay class inventories to pixel-level biomass (Mg ha⁻¹) using density reduction factors |

## Usage

```r
# Install SpaDES.project if needed
install.packages("SpaDES.project", repos = "https://predictiveecology.r-universe.dev")

source("global.R")
```

`global.R` will:

1. Install any missing R packages via `Require`.
2. Download the three modules from GitHub via `SpaDES.install::installModules()`.
3. Run a 50-year example simulation on a 3×3 pixel *Pinus strobus* landscape.
4. Print a biomass summary and display spatial and time-series plots.

## Repository layout

```
SpaDES_PinStro_DeadWoodExample/
├── global.R            # entry point — run this
├── R/
│   ├── parameters.R    # transition matrices, fall/loss probabilities, DRF lookup
│   └── example-data.R  # cohort data and study area raster for the example run
├── modules/            # populated at runtime by installModules() — gitignored
├── inputs/             # simulation inputs — gitignored
├── outputs/            # simulation outputs — gitignored
└── cache/              # reproducible::Cache() store — gitignored
```
