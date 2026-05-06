library(SpaDES.core)
library(data.table)
library(terra)
library(testthat)

DRFLookup_test <- data.table(
  species = "Pinus strobus",
  pool    = rep(c("snag", "DWD"), each = 5),
  DC      = rep(1:5, times = 2),
  DRF     = c(1.000, 0.841, 0.706, 0.543, 0.382,
              1.000, 0.783, 0.614, 0.418, 0.251)
)

templateRaster <- terra::rast(nrows = 3, ncols = 3,
                               xmin = 0, xmax = 3, ymin = 0, ymax = 3,
                               crs = "EPSG:4326")

emptySnagTable <- data.table(
  pixelID = integer(), species = character(),
  DC = integer(), ageInDC = integer(), initBiomass = numeric()
)

# testInit is not available in SpaDES.core >= 3.x; use simInit directly.
# This test file lives at: modules/deadWoodBiomass/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}

testInit <- function(moduleName, params, objects, times = list(start = 0, end = 10)) {
  simInit(
    times   = times,
    modules = list(moduleName),
    params  = params,
    objects = objects,
    paths   = list(modulePath = file.path(.projRoot, "modules"))
  )
}

test_that("deadWoodBiomass init creates snagBiomass and DWDBiomass rasters", {
  sim <- testInit(
    "deadWoodBiomass",
    params  = list(deadWoodBiomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = data.table::copy(emptySnagTable),
      DWDTable        = data.table::copy(emptySnagTable),
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = "init")
  expect_s4_class(sim$snagBiomass_Mg_ha, "SpatRaster")
  expect_s4_class(sim$DWDBiomass_Mg_ha,  "SpatRaster")
  expect_true(all(is.na(terra::values(sim$snagBiomass_Mg_ha))))
  expect_true(all(is.na(terra::values(sim$DWDBiomass_Mg_ha))))
})
