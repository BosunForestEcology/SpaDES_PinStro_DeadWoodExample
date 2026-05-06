library(SpaDES.project)  # bootstrap exception — only library() call permitted

Require::Require(c("SpaDES.core", "data.table", "terra"))

source("R/parameters.R")

# Minimal single-pixel cohort data for a 50-year run
myMortalityTable <- data.table::data.table(
  pixelID = 1L,
  year    = c(1L, 5L, 10L),    # mortality events at years 1, 5, 10
  species = "Pinus strobus",
  B       = c(15.0, 10.0, 8.0) # Mg ha-1
)

myRaster <- terra::rast(
  nrows = 1, ncols = 1,
  xmin = 0, xmax = 1, ymin = 0, ymax = 1,
  crs = "EPSG:4326"  # geometry-only template; values are set by deadWoodBiomass
)

times   <- list(start = 0, end = 50)
params  <- list(
  snagDecay = list(
    snagTransMat = snagTransMat,
    snagFallProb = snagFallProb,
    species      = "Pinus strobus"
  ),
  DWDDecay = list(
    DWDTransMat     = DWDTransMat,
    snagToDWD_DCmap = snagToDWD_DCmap,
    DWD_lossProb    = DWD_lossProb
  ),
  deadWoodBiomass = list(
    DRFLookup        = DRFLookup,
    .plotInitialTime = NA
  )
)
modules <- list("snagDecay", "DWDDecay", "deadWoodBiomass")

# Note: fallenSnags is not provided here; snagDecayInit() creates it at time 0
# (before DWDDecay's first receive event at time 1), so the SpaDES contract
# warning for fallenSnags at simInit is expected and resolves correctly at runtime.
mySim <- simInit(
  times   = times,
  params  = params,
  modules = modules,
  objects = list(
    cohortData      = myMortalityTable,
    studyAreaRaster = myRaster
  ),
  paths = list(
    modulePath = file.path("modules"),
    inputPath  = file.path("inputs"),
    outputPath = file.path("outputs"),
    cachePath  = file.path("cache")
  )
)

set.seed(42)
mySim <- spades(mySim)

# Inspect outputs
cat("Final snag inventory rows:   ", nrow(mySim$snagTable), "\n")
cat("Final DWD inventory rows:    ", nrow(mySim$DWDTable),  "\n")
cat("Final snag biomass (pixel 1):", terra::values(mySim$snagBiomass_Mg_ha)[1, 1], "Mg/ha\n")
cat("Final DWD biomass  (pixel 1):", terra::values(mySim$DWDBiomass_Mg_ha)[1, 1],  "Mg/ha\n")
