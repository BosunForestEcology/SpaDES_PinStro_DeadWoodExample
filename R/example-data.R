# R/example-data.R
# Generates a 3x3 (9-pixel) study area with staggered mortality events
# for a 50-year Pinus strobus dead wood decay run.
# Source this file from global.R before calling simInit().

# ---- Cohort data --------------------------------------------------------
# Mortality enters each pixel at different years, with varying biomass loads
# to exercise the full pipeline across pixels.
myMortalityTable <- data.table::data.table(
  pixelID = c(
    1L, 1L, 1L,
    2L, 2L,
    3L, 3L,
    4L,
    5L, 5L, 5L,
    6L, 6L,
    7L,
    8L, 8L,
    9L
  ),
  year = c(
    1L,  8L, 20L,   # pixel 1: early, mid, late pulses
    3L, 15L,        # pixel 2
    5L, 25L,        # pixel 3
    2L,             # pixel 4: single early pulse
    1L, 10L, 30L,   # pixel 5: spread over run
    7L, 18L,        # pixel 6
   12L,             # pixel 7: mid-run only
    4L, 22L,        # pixel 8
   35L              # pixel 9: late entry
  ),
  species     = "Pinus strobus",
  diameter_cm = c(
    18.2, 14.5,  9.8,
    22.1, 16.3,
    20.4, 11.7,
    25.0,
    13.6, 17.8,  8.9,
    24.3, 15.1,
    19.5,
    16.2, 10.4,
    12.8
  ),
  biomass = c(
    20.0, 12.0,  8.0,
    15.0, 10.0,
    18.0,  9.0,
    25.0,
    10.0, 14.0,  6.0,
    22.0, 11.0,
    16.0,
    13.0,  7.0,
     9.0
  )
)

# ---- Study area raster --------------------------------------------------
# 3 x 3 grid, 100 m resolution, arbitrary location in EPSG:32610 (UTM 10N)
myRaster <- terra::rast(
  nrows = 3, ncols = 3,
  xmin = 0, xmax = 300,
  ymin = 0, ymax = 300,
  crs = "EPSG:32610",
  resolution = 100
)
