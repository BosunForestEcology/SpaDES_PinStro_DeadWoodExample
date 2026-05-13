# R/example-data.R
# Generates a 6x3 (18-pixel) study area with staggered mortality events
# for a 100-year Pinus strobus + Pinus resinosa dead wood decay run.
# Pixels 1-9 occupy the top half of the grid (same positions as the original
# 3x3 layout); pixels 10-18 are the new bottom half.
# All mortality events occur at year <= 30.
# Source this file from global.R before calling simInit().

# ---- Cohort data --------------------------------------------------------
# Pixels 1-6: Pinus strobus; Pixels 4-9: Pinus resinosa (pixels 4-6 have both).
# Pixels 10-18: mixed species, staggered years, varied diameters and biomass.
myMortalityTable <- data.table::data.table(
  pixelID = c(
    # Pinus strobus — original pixels
    1L, 1L, 1L,
    2L, 2L,
    3L, 3L,
    4L,
    5L, 5L,
    6L, 6L,
    # Pinus strobus — new pixels
    10L,
    11L, 11L,
    13L,
    14L,
    16L,
    17L,
    18L,
    # Pinus resinosa — original pixels
    4L,
    5L,
    6L, 6L,
    7L, 7L,
    8L, 8L,
    9L,
    # Pinus resinosa — new pixels
    10L,
    12L, 12L,
    14L,
    15L, 15L,
    16L,
    17L
  ),
  year = c(
    # Pinus strobus — original pixels
    1L,  8L, 20L,
    3L, 15L,
    5L, 25L,
    2L,
    1L, 30L,
    7L, 18L,
    # Pinus strobus — new pixels
    5L,
    12L, 28L,
    9L,
    2L,
    16L,
    23L,
    11L,
    # Pinus resinosa — original pixels
    10L,
    12L,
     4L, 22L,
     6L, 28L,
     3L, 28L,
    15L,
    # Pinus resinosa — new pixels
    18L,
     4L, 22L,
    20L,
     7L, 30L,
     6L,
    14L
  ),
  species = c(
    rep("Pinus strobus",  20L),
    rep("Pinus resinosa", 17L)
  ),
  diameter_cm = c(
    # Pinus strobus — original pixels
    18.2, 14.5,  9.8,
    22.1, 16.3,
    20.4, 11.7,
    25.0,
    13.6, 17.8,
    24.3, 15.1,
    # Pinus strobus — new pixels
    19.0,
    16.5, 10.5,
    21.0,
    24.5,
    18.6,
    13.8,
    20.5,
    # Pinus resinosa — original pixels
    21.5,
    19.2,
    23.0, 12.4,
    17.8, 10.9,
    20.1, 14.6,
    16.3,
    # Pinus resinosa — new pixels
    22.4,
    23.8, 14.2,
    15.6,
    20.0, 11.8,
    21.3,
    17.0
  ),
  biomass = c(
    # Pinus strobus — original pixels
    20.0, 12.0,  8.0,
    15.0, 10.0,
    18.0,  9.0,
    25.0,
    10.0, 17.0,
    22.0, 11.0,
    # Pinus strobus — new pixels
    16.0,
    13.0,  8.0,
    18.0,
    23.0,
    14.0,
    11.0,
    19.0,
    # Pinus resinosa — original pixels
    18.0,
    14.0,
    20.0,  8.0,
    16.0,  7.0,
    19.0, 11.0,
    13.0,
    # Pinus resinosa — new pixels
    15.0,
    20.0, 10.0,
    12.0,
    17.0,  7.0,
    18.0,
    13.0
  )
)

# ---- Study area raster --------------------------------------------------
# 6x3 grid, 100 m resolution, arbitrary location in EPSG:32610 (UTM 10N).
# Pixels 1-9 (top half) match the original 3x3 layout.
myRaster <- terra::rast(
  nrows = 6, ncols = 3,
  xmin = 0, xmax = 300,
  ymin = 0, ymax = 600,
  crs = "EPSG:32610",
  resolution = 100
)
