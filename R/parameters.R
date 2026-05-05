# R/parameters.R
# All shared decay parameters for the White Pine dead wood model.
# Source this file from global.R before calling simInit().

snagTransMat <- matrix(
  c(
    0.48, 0.38, 0.05, 0.00, 0.00,  # from DC1
    0.00, 0.52, 0.34, 0.05, 0.00,  # from DC2
    0.00, 0.00, 0.56, 0.30, 0.04,  # from DC3
    0.00, 0.00, 0.00, 0.62, 0.26,  # from DC4
    0.00, 0.00, 0.00, 0.00, 0.70   # from DC5
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))
)
# NOTE: Replace placeholder values above with exact values from Paper 1, Table 2.

snagFallProb <- c(DC1 = 0.09, DC2 = 0.09, DC3 = 0.10, DC4 = 0.12, DC5 = 0.30)
# NOTE: Replace with exact species-specific values from Paper 1.

DWDTransMat <- matrix(
  c(
    0.55, 0.35, 0.06, 0.00, 0.00,  # from DC1
    0.00, 0.50, 0.37, 0.08, 0.00,  # from DC2
    0.00, 0.00, 0.48, 0.38, 0.08,  # from DC3
    0.00, 0.00, 0.00, 0.50, 0.36,  # from DC4
    0.00, 0.00, 0.00, 0.00, 0.72   # from DC5
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))
)
# NOTE: Replace placeholder values above with exact values from Paper 1, Table 2.
# NOTE: DC1 row sums to 0.96 with DWD_lossProb[DC1] = 0 — verify exact value against paper.

DWD_lossProb <- c(DC1 = 0.00, DC2 = 0.05, DC3 = 0.06, DC4 = 0.14, DC5 = 0.28)
# NOTE: Replace with exact species-specific values from Paper 1.

snagToDWD_DCmap <- c(DC1 = 1L, DC2 = 2L, DC3 = 2L, DC4 = 3L, DC5 = 4L)
# NOTE: Confirm DC mapping against Paper 1 methods. DC5 snags enter DWD at DC4 (not DC5).

DRFLookup <- data.table::data.table(
  species = "Pinus strobus",
  pool    = rep(c("snag", "DWD"), each = 5),
  DC      = rep(1:5, times = 2),
  DRF     = c(
    1.000, 0.841, 0.706, 0.543, 0.382,  # snag
    1.000, 0.783, 0.614, 0.418, 0.251   # DWD
  )
)
# NOTE: Replace placeholder DRF values above with exact values from Paper 2, Appendix D.
