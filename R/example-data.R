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

# ---- Pre-existing snag inventory ----------------------------------------
# ~50 rows of standing dead trees already present at simulation start.
# DC1-DC5 are all represented; the proportion of higher DCs increases toward
# the lower rows (older, more established forest).
snag_pixels <- c(
   5L,  8L, 12L, 15L, 18L, 22L, 25L, 30L, 33L, 37L,   # row 1-2 (young)
  42L, 48L, 55L, 60L, 65L, 70L, 75L, 82L, 88L, 95L,   # row 3-5 (mixed)
  101L, 108L, 115L, 122L, 130L, 138L, 145L, 153L,       # row 6-8 (older)
  161L, 168L, 175L, 182L, 188L, 195L, 198L              # row 9-10 (established)
)
snag_species <- rep(c("Pinus strobus", "Pinus resinosa"),
                    length.out = length(snag_pixels))

row_of_snag <- ((snag_pixels - 1L) %/% 20L) + 1L

dc_sample <- vapply(row_of_snag, function(r) {
  p_dc <- c(0.30, 0.25, 0.20, 0.15, 0.10) +
          c(-0.04, -0.02, 0.01, 0.03, 0.02) * (r - 1L)
  p_dc <- pmax(p_dc, 0.01)
  p_dc <- p_dc / sum(p_dc)
  sample(1:5, 1L, prob = p_dc)
}, integer(1L))

myInitialSnagData <- data.table::data.table(
  pixelID     = snag_pixels,
  species     = snag_species,
  DC          = as.integer(dc_sample),
  ageInDC     = as.integer(round(stats::runif(length(snag_pixels), 0, 20))),
  initBiomass = round(stats::runif(length(snag_pixels), 18, 45), 1),
  diameter_cm = round(stats::runif(length(snag_pixels), 20, 40), 1)
)
data.table::setorder(myInitialSnagData, pixelID, species)

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
