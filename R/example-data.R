# R/example-data.R
# Generates a 3x3 (9-pixel) study area with staggered mortality events
# for a 100-year Pinus strobus + Pinus resinosa dead wood decay run.
# Source this file from global.R before calling simInit().

# ---- Cohort data --------------------------------------------------------
# Pixels 1-6: Pinus strobus; Pixels 4-9: Pinus resinosa (pixels 4-6 have both).
# Mortality enters at different years and biomass loads to exercise the full pipeline.
myMortalityTable <- data.table::data.table(
  pixelID = c(
    # Pinus strobus
    1L, 1L, 1L,
    2L, 2L,
    3L, 3L,
    4L,
    5L, 5L,
    6L, 6L,
    # Pinus resinosa
    4L,
    5L,
    6L, 6L,
    7L, 7L,
    8L, 8L,
    9L
  ),
  year = c(
    # Pinus strobus
    1L,  8L, 20L,
    3L, 15L,
    5L, 25L,
    2L,
    1L, 30L,
    7L, 18L,
    # Pinus resinosa
    10L,
    12L,
     4L, 22L,
     6L, 28L,
     3L, 35L,
    15L
  ),
  species = c(
    rep("Pinus strobus",  12L),
    rep("Pinus resinosa",  9L)
  ),
  diameter_cm = c(
    # Pinus strobus
    18.2, 14.5,  9.8,
    22.1, 16.3,
    20.4, 11.7,
    25.0,
    13.6, 17.8,
    24.3, 15.1,
    # Pinus resinosa
    21.5,
    19.2,
    23.0, 12.4,
    17.8, 10.9,
    20.1, 14.6,
    16.3
  ),
  biomass = c(
    # Pinus strobus
    20.0, 12.0,  8.0,
    15.0, 10.0,
    18.0,  9.0,
    25.0,
    10.0, 17.0,
    22.0, 11.0,
    # Pinus resinosa
    18.0,
    14.0,
    20.0,  8.0,
    16.0,  7.0,
    19.0, 11.0,
    13.0
  )
)

# ---- Study area raster --------------------------------------------------
# 3x3 grid, 100 m resolution, arbitrary location in EPSG:32610 (UTM 10N)
myRaster <- terra::rast(
  nrows = 3, ncols = 3,
  xmin = 0, xmax = 300,
  ymin = 0, ymax = 300,
  crs = "EPSG:32610",
  resolution = 100
)
