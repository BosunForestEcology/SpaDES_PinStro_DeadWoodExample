# R/example-data.R
# Generates a 20×10 (200-pixel) live forest cohort for a 200-year dead wood decay run.
# Both species may occupy the same pixel. Age-class presence follows a spatial gradient:
# young cohorts are more common in upper rows (recent disturbance origin) and old cohorts
# are more common in lower rows (established forest). Pixel IDs run left-to-right,
# top-to-bottom (terra default), matching the 20-column raster below.

set.seed(42)

n_pixels <- 200L
species  <- c("Pinus strobus", "Pinus resinosa")

# Row index for each pixel (row 1 = top, row 10 = bottom)
pixel_row <- ((seq_len(n_pixels) - 1L) %/% 20L) + 1L

rows <- lapply(seq_len(n_pixels), function(px) {
  row_idx <- pixel_row[px]

  # Age-class presence probabilities shift with row
  p_young <- 0.70 - 0.04 * row_idx  # more young trees near disturbance origin (top)
  p_mid   <- 0.50
  p_old   <- 0.20 + 0.05 * row_idx  # more old trees in established forest (bottom)

  dt_list <- list()
  for (sp in species) {
    if (stats::runif(1L) < p_young) {
      dt_list[[length(dt_list) + 1L]] <- data.table::data.table(
        pixelID     = px,
        species     = sp,
        biomass     = round(stats::runif(1L, 8,  15), 1),
        diameter_cm = round(stats::runif(1L, 10, 18), 1)
      )
    }
    if (stats::runif(1L) < p_mid) {
      dt_list[[length(dt_list) + 1L]] <- data.table::data.table(
        pixelID     = px,
        species     = sp,
        biomass     = round(stats::runif(1L, 18, 28), 1),
        diameter_cm = round(stats::runif(1L, 20, 28), 1)
      )
    }
    if (stats::runif(1L) < p_old) {
      dt_list[[length(dt_list) + 1L]] <- data.table::data.table(
        pixelID     = px,
        species     = sp,
        biomass     = round(stats::runif(1L, 30, 45), 1),
        diameter_cm = round(stats::runif(1L, 30, 40), 1)
      )
    }
  }
  if (length(dt_list) == 0L) return(NULL)
  data.table::rbindlist(dt_list)
})

myLiveCohortData <- data.table::rbindlist(Filter(Negate(is.null), rows))
data.table::setorder(myLiveCohortData, pixelID, species)

# ---- Study area raster --------------------------------------------------
# 20 cols × 10 rows, 10 m resolution — 200 pixels total.
myRaster <- terra::rast(
  nrows      = 10,
  ncols      = 20,
  xmin       = 0,
  xmax       = 200,
  ymin       = 0,
  ymax       = 100,
  crs        = "EPSG:32610",
  resolution = 10
)
