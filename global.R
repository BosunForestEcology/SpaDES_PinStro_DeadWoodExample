repos <- unique(c("predictiveecology.r-universe.dev", getOption("repos")))
install.packages("SpaDES.project", repos = repos)

projectPath <- getwd()

source("R/example-data.R")  # myLiveCohortData, myRaster

times <- list(start = 0, end = 200)

# Note: fallenSnags is not provided here; DeadWood_snagDecay Init() creates it
# at time 0 (before DeadWood_DWDDecay's first receive event at time 5), so the
# SpaDES contract warning for fallenSnags at simInit is expected and resolves
# correctly at runtime.
out <- SpaDES.project::setupProject(
  useGit  = FALSE,
  overwrite = TRUE,
  paths   = list(
    projectPath = projectPath,
    modulePath  = file.path(projectPath, "modules"),
    inputPath   = file.path(projectPath, "inputs"),
    outputPath  = file.path(projectPath, "outputs"),
    cachePath   = file.path(projectPath, "cache")
  ),
  options = options(
    repos                        = c(repos = repos),
    reproducible.destinationPath = "inputs",
    reproducible.useMemoise      = TRUE,
    spades.moduleCodeChecks      = FALSE
  ),
  modules = c(
    "BosunForestEcology/DeadWood_Mortality@main",
    "BosunForestEcology/DeadWood_snagDecay@main",
    "BosunForestEcology/DeadWood_DWDDecay@main",
    "BosunForestEcology/DeadWood_Biomass@main"
  ),
  times  = times,
  params = list(
    DeadWood_Mortality = list(baseMortality = 0.001),
    DeadWood_snagDecay = list(species = c("Pinus strobus", "Pinus resinosa")),
    DeadWood_Biomass   = list(.plotInitialTime = 5)
  ),
  liveCohortData  = myLiveCohortData,
  initialSnagTable = myInitialSnagData,
  studyAreaRaster = myRaster
)

#set.seed(42)
mySim <- SpaDES.core::simInitAndSpades2(out)

# Inspect outputs
cat("Final snag inventory rows:   ", nrow(mySim$snagTable), "\n")
cat("Final DWD inventory rows:    ", nrow(mySim$DWDTable),  "\n")
cat("Final snag biomass (pixel 1):", terra::values(mySim$snagBiomass_Mg_ha)[1, 1], "Mg/ha\n")
cat("Final DWD biomass  (pixel 1):", terra::values(mySim$DWDBiomass_Mg_ha)[1, 1],  "Mg/ha\n")

# ---- Visualize biomass history ------------------------------------------

# 1. Spatial snapshots — per species, grid at each 5-year step
for (sp in names(mySim$snagHistoryBySpecies)) {
  snap_years <- sub("yr", "Year ", names(mySim$snagHistoryBySpecies[[sp]]))
  snag_range <- range(terra::values(mySim$snagHistoryBySpecies[[sp]]), na.rm = TRUE)
  DWD_range  <- range(terra::values(mySim$DWDHistoryBySpecies[[sp]]),  na.rm = TRUE)

  terra::plot(mySim$snagHistoryBySpecies[[sp]],
              main   = paste0(sp, " — Snag — ", snap_years),
              range  = snag_range,
              col    = hcl.colors(50, "YlOrRd", rev = TRUE),
              legend = "bottomright")

  terra::plot(mySim$DWDHistoryBySpecies[[sp]],
              main   = paste0(sp, " — DWD — ", snap_years),
              range  = DWD_range,
              col    = hcl.colors(50, "Greens", rev = TRUE),
              legend = "bottomright")
}

# 2. Time-series line chart — total biomass per pool per species, faceted by species
bioSummary <- data.table::rbindlist(lapply(names(mySim$snagHistoryBySpecies), function(sp) {
  hist_s <- mySim$snagHistoryBySpecies[[sp]]
  hist_d <- mySim$DWDHistoryBySpecies[[sp]]
  years  <- as.integer(sub("yr", "", names(hist_s)))
  data.table::data.table(
    year    = rep(years, 2L),
    pool    = rep(c("Snag", "DWD"), each = length(years)),
    species = sp,
    total   = c(colSums(terra::values(hist_s), na.rm = TRUE),
                colSums(terra::values(hist_d),  na.rm = TRUE))
  )
}))

print(
  ggplot2::ggplot(bioSummary, ggplot2::aes(x = year, y = total, colour = pool)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(DWD = "#2c8c4f", Snag = "#c0392b")) +
    ggplot2::facet_wrap(~ species) +
    ggplot2::labs(
      title  = "Dead wood biomass over time (200 pixels)",
      x      = "Year",
      y      = "Total biomass (Mg/ha, summed across pixels)",
      colour = "Pool"
    ) +
    ggplot2::theme_bw(base_size = 13)
)
