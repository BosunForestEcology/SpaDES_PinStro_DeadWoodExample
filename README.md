# DeadWood_Example

Example repository for a *Pinus strobus* and *Pinus resinosa* snag and dead wood decay workflow built
with [SpaDES](https://spades.predictiveecology.org/).

These modules are built on top of the work from Vanderwel, M.C., Malcolm, J.R., Smith, S.M., and Islam, N. (2006). An
integrated model for snag and downed woody debris decay class
transition. *Forest Ecology and Management*, 234(1–3), 48–59.
<https://doi.org/10.1016/j.foreco.2006.06.020>

Clone this repository and then run `global.R`. 

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
3. Run a 100-year example simulation on a 3×3 pixel *Pinus strobus* and *Pinus resinosa* landscape.
4. Print a biomass summary and display spatial and time-series plots.

## Repository layout

```
DeadWood_Example/
├── global.R            # entry point — run this
├── R/
│   ├── parameters.R    # transition matrices, fall/loss probabilities, DRF lookup
│   └── example-data.R  # cohort data and study area raster for the example run
├── modules/            # populated at runtime by installModules() — gitignored
├── inputs/             # simulation inputs — gitignored
├── outputs/            # simulation outputs — gitignored
└── cache/              # reproducible::Cache() store — gitignored
```
