# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "deadWoodBiomass",
  description = "Translates snag and DWD decay class inventories into pixel-level biomass
                 estimates (Mg ha-1) using species- and pool-specific density reduction
                 factors (DRF) from Paper 2 Appendix D.",
  keywords    = c("dead wood", "biomass", "density reduction factor", "carbon"),
  authors     = structure(list(list(given = "First", family = "Last",
                                    role = c("aut", "cre"),
                                    email = "email@example.com", comment = NULL)),
                           class = "person"),
  childModules = character(0),
  version     = list(deadWoodBiomass = "0.0.1"),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "terra", "SpaDES.core (>= 3.0.0)"),
  parameters  = bindrows(
    defineParameter("DRFLookup", "data.table",
                    data.table::data.table(species = character(), pool = character(),
                                           DC = integer(), DRF = numeric()),
                    NA, NA,
                    desc = "Density reduction factor lookup: species, pool, DC, DRF."),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "Simulation time for first plot. NA = no plots."),
    defineParameter(".plotInterval",    "numeric",  1, NA, NA,
                    desc = "Interval between plots.")
  ),
  inputObjects = bindrows(
    expectsInput("snagTable", "data.table",
                 desc = "Current snag inventory from snagDecay."),
    expectsInput("DWDTable", "data.table",
                 desc = "Current DWD inventory from DWDDecay."),
    expectsInput("studyAreaRaster", "SpatRaster",
                 desc = "Template raster defining pixel grid, CRS, and resolution.")
  ),
  outputObjects = bindrows(
    createsOutput("snagBiomass_Mg_ha", "SpatRaster",
                  desc = "Pixel-level snag biomass (Mg ha-1)."),
    createsOutput("DWDBiomass_Mg_ha", "SpatRaster",
                  desc = "Pixel-level DWD biomass (Mg ha-1).")
  )
))

doEvent.deadWoodBiomass <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- deadWoodBiomassInit(sim)
      sim <- scheduleEvent(sim, start(sim) + 1, "deadWoodBiomass", "annual", eventPriority = 4)
    },
    annual = {
      sim <- deadWoodBiomassAnnual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "deadWoodBiomass", "annual", eventPriority = 4)
    },
    warning(paste("Undefined event type:", eventType, "in module deadWoodBiomass"))
  )
  return(invisible(sim))
}

deadWoodBiomassInit <- function(sim) {
  if (nrow(P(sim)$DRFLookup) == 0L)
    stop("DRFLookup is empty — provide a density reduction factor table in params.")
  # terra::values<- deep-copies before writing, so studyAreaRaster is not modified
  sim$snagBiomass_Mg_ha <- sim$studyAreaRaster
  terra::values(sim$snagBiomass_Mg_ha) <- NA_real_
  sim$DWDBiomass_Mg_ha  <- sim$studyAreaRaster
  terra::values(sim$DWDBiomass_Mg_ha)  <- NA_real_
  return(invisible(sim))
}

deadWoodBiomassAnnual <- function(sim) {
  drf <- P(sim)$DRFLookup

  snagWithBiomass <- data.table::copy(sim$snagTable)
  snagWithBiomass <- drf[pool == "snag"][snagWithBiomass, on = .(species, DC)]
  snagWithBiomass[, currentBiomass := initBiomass * DRF]
  snagByPixel <- snagWithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
  sim$snagBiomass_Mg_ha <- pixelValuesToRaster(snagByPixel, sim$studyAreaRaster)

  DWDwithBiomass <- data.table::copy(sim$DWDTable)
  DWDwithBiomass <- drf[pool == "DWD"][DWDwithBiomass, on = .(species, DC)]
  DWDwithBiomass[, currentBiomass := initBiomass * DRF]
  DWDbyPixel <- DWDwithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
  sim$DWDBiomass_Mg_ha <- pixelValuesToRaster(DWDbyPixel, sim$studyAreaRaster)

  return(invisible(sim))
}
