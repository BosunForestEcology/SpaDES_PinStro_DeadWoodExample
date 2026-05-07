repos <- unique(c("predictiveecology.r-universe.dev", getOption("repos")))
install.packages("SpaDES.project", repos = repos)
Require::Require(c("SpaDES.core", "data.table", "terra", "ggplot2"))

projectPath <- getwd()
modulePath  <- file.path(projectPath, "modules")

for (d in c(modulePath, file.path(projectPath, "inputs"),
            file.path(projectPath, "outputs"), file.path(projectPath, "cache"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Download modules from BosunForestEcology GitHub as zip archives
for (mod in c("DeadWood_snagDecay", "DeadWood_DWDDecay", "DeadWood_Biomass")) {
  modPath <- file.path(modulePath, mod)
  if (!dir.exists(modPath)) {
    tmp <- tempfile(fileext = ".zip")
    download.file(
      url      = paste0("https://github.com/BosunForestEcology/", mod, "/archive/refs/heads/main.zip"),
      destfile = tmp,
      mode     = "wb",
      quiet    = TRUE
    )
    unzip(tmp, exdir = modulePath)
    file.rename(file.path(modulePath, paste0(mod, "-main")), modPath)
    unlink(tmp)
  }
}

source("R/example-data.R")  # myMortalityTable, myRaster

times <- list(start = 0, end = 50)

# Note: fallenSnags is not provided here; DeadWood_snagDecay Init() creates it
# at time 0 (before DeadWood_DWDDecay's first receive event at time 5), so the
# SpaDES contract warning for fallenSnags at simInit is expected and resolves
# correctly at runtime.
mySim <- SpaDES.core::simInit(
  times   = times,
  params  = list(DeadWood_Biomass = list(.plotInitialTime = 5)),
  modules = list("DeadWood_snagDecay", "DeadWood_DWDDecay", "DeadWood_Biomass"),
  objects = list(cohortData = myMortalityTable, studyAreaRaster = myRaster),
  paths   = list(
    modulePath = modulePath,
    inputPath  = file.path(projectPath, "inputs"),
    outputPath = file.path(projectPath, "outputs"),
    cachePath  = file.path(projectPath, "cache")
  )
)

set.seed(42)
mySim <- SpaDES.core::spades(mySim)

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
