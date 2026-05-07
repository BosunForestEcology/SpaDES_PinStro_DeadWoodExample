library(SpaDES.project)  # bootstrap exception — only library() call permitted

Require::Require(c("SpaDES.core", "SpaDES.install", "data.table", "terra", "ggplot2"))

source("R/parameters.R")
source("R/example-data.R")   # myMortalityTable, myRaster

for (d in c("inputs", "outputs", "cache")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Download modules from BosunForestEcology if not already present
SpaDES.install::installModules(
  modules    = c("DeadWood_snagDecay", "DeadWood_DWDDecay", "DeadWood_Biomass"),
  modulePath = "modules",
  account    = "BosunForestEcology"
)

times   <- list(start = 0, end = 50)
params  <- list(
  DeadWood_snagDecay = list(
    snagTransMat = snagTransMat,
    snagFallProb = snagFallProb,
    species      = "Pinus strobus"
  ),
  DeadWood_DWDDecay = list(
    DWDTransMat     = DWDTransMat,
    snagToDWD_DCmap = snagToDWD_DCmap,
    DWD_lossProb    = DWD_lossProb
  ),
  DeadWood_Biomass = list(
    DRFLookup        = DRFLookup,
    .plotInitialTime = 5,
    .plotInterval    = 5
  )
)
modules <- list("DeadWood_snagDecay", "DeadWood_DWDDecay", "DeadWood_Biomass")

# Note: fallenSnags is not provided here; DeadWood_snagDecay Init() creates it
# at time 0 (before DeadWood_DWDDecay's first receive event at time 1), so the
# SpaDES contract warning for fallenSnags at simInit is expected and resolves
# correctly at runtime.
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

# ---- Visualize biomass history ------------------------------------------

# 1. Spatial snapshots — 3x3 grid at each 5-year step
snag_range <- range(terra::values(mySim$snagHistory), na.rm = TRUE)
DWD_range  <- range(terra::values(mySim$DWDHistory),  na.rm = TRUE)
snap_years <- sub("yr", "Year ", names(mySim$snagHistory))

terra::plot(mySim$snagHistory,
            main   = paste("Snag biomass —", snap_years),
            range  = snag_range,
            col    = hcl.colors(50, "YlOrRd", rev = TRUE),
            legend = "bottomright")

terra::plot(mySim$DWDHistory,
            main   = paste("DWD biomass —", snap_years),
            range  = DWD_range,
            col    = hcl.colors(50, "Greens", rev = TRUE),
            legend = "bottomright")

# 2. Time-series line chart — total biomass per pool across all pixels
years      <- as.integer(sub("yr", "", names(mySim$snagHistory)))
snagTotals <- colSums(terra::values(mySim$snagHistory), na.rm = TRUE)
DWDTotals  <- colSums(terra::values(mySim$DWDHistory),  na.rm = TRUE)

bioSummary <- data.table::data.table(
  year  = rep(years, 2L),
  pool  = rep(c("Snag", "DWD"), each = length(years)),
  total = c(snagTotals, DWDTotals)
)

print(
  ggplot2::ggplot(bioSummary, ggplot2::aes(x = year, y = total, colour = pool)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(DWD = "#2c8c4f", Snag = "#c0392b")) +
    ggplot2::labs(
      title  = "Dead wood biomass over time (9 pixels, Pinus strobus)",
      x      = "Year",
      y      = "Total biomass (Mg/ha, summed across pixels)",
      colour = "Pool"
    ) +
    ggplot2::theme_bw(base_size = 13)
)
