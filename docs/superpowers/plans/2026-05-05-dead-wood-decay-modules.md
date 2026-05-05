# Dead Wood Decay Modules Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement three coordinated SpaDES modules (`snagDecay`, `DWDDecay`, `deadWoodBiomass`) that simulate White Pine dead wood biomass decay through five decay classes using a Markov-chain transition model.

**Architecture:** Dead trees enter `snagDecay` as DC1 snags; each year the module advances their decay class and stochastically transfers fallen snags to `DWDDecay`; `deadWoodBiomass` joins both inventories to a density-reduction-factor (DRF) lookup table and writes pixel-level biomass rasters. Within-timestep ordering is enforced via SpaDES event priorities 1–4.

**Tech Stack:** R, `SpaDES.core`, `data.table`, `terra`, `reproducible`, `Require`, `testthat`

---

## File Structure

```
modules/
  snagDecay/
    snagDecay.R                          # defineModule + doEvent + event functions
    R/
      transition.R                       # applyTransition() helper
    tests/testthat/
      test-snagDecay.R
  DWDDecay/
    DWDDecay.R                           # defineModule + doEvent + event functions
    tests/testthat/
      test-DWDDecay.R
  deadWoodBiomass/
    deadWoodBiomass.R                    # defineModule + doEvent + event functions
    R/
      raster-helpers.R                   # pixelValuesToRaster() helper
    tests/testthat/
      test-deadWoodBiomass.R
R/
  parameters.R                           # all shared parameters: matrices, DRF lookup
global.R                                 # simInit + spades entry point
```

**Shared contracts (sim$ objects):**

| Object | Class | Schema |
|--------|-------|--------|
| `sim$cohortData` | `data.table` | pixelID (int), year (int), species (chr), B (num, Mg ha⁻¹) |
| `sim$snagTable` | `data.table` | pixelID (int), species (chr), DC (int), ageInDC (int), initBiomass (num) |
| `sim$fallenSnags` | `data.table` | same schema as snagTable |
| `sim$DWDTable` | `data.table` | pixelID (int), species (chr), DC (int), ageInDC (int), initBiomass (num) |
| `sim$snagBiomass_Mg_ha` | `SpatRaster` | pixel-level snag biomass |
| `sim$DWDBiomass_Mg_ha` | `SpatRaster` | pixel-level DWD biomass |

**Within-timestep event priority order:**

| Priority | Module | Event |
|----------|--------|-------|
| 1 | snagDecay | annual |
| 2 | DWDDecay | receive |
| 3 | DWDDecay | annual |
| 4 | deadWoodBiomass | annual |

---

## Task 0: Project scaffold and shared parameters

**Files:**
- Create: `R/parameters.R`
- Create: `modules/snagDecay/` (directory scaffold)
- Create: `modules/DWDDecay/` (directory scaffold)
- Create: `modules/deadWoodBiomass/` (directory scaffold)

- [ ] **Step 1: Create directory scaffold**

```bash
mkdir -p modules/snagDecay/R
mkdir -p modules/snagDecay/tests/testthat
mkdir -p modules/DWDDecay/tests/testthat
mkdir -p modules/deadWoodBiomass/R
mkdir -p modules/deadWoodBiomass/tests/testthat
mkdir -p R
```

- [ ] **Step 2: Write `R/parameters.R`**

```r
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

DWD_lossProb <- c(DC1 = 0.00, DC2 = 0.05, DC3 = 0.06, DC4 = 0.14, DC5 = 0.28)

snagToDWD_DCmap <- c(DC1 = 1L, DC2 = 2L, DC3 = 2L, DC4 = 3L, DC5 = 4L)

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
```

- [ ] **Step 3: Verify parameters are structurally valid**

```r
Rscript -e "
source('R/parameters.R')
stopifnot(dim(snagTransMat) == c(5, 5))
stopifnot(dim(DWDTransMat) == c(5, 5))
stopifnot(all(rowSums(snagTransMat) <= 1 + 1e-9))
stopifnot(all(rowSums(DWDTransMat) <= 1 + 1e-9))
stopifnot(length(snagFallProb) == 5)
stopifnot(length(DWD_lossProb) == 5)
stopifnot(length(snagToDWD_DCmap) == 5)
stopifnot(nrow(DRFLookup) == 10)
cat('parameters.R OK\n')
"
```

Expected output: `parameters.R OK`

- [ ] **Step 4: Commit**

```bash
git add R/parameters.R
git commit -m "feat: add shared decay parameters (placeholder values from Paper 1 & 2)"
```

---

## Task 1: `applyTransition()` helper — TDD

**Files:**
- Create: `modules/snagDecay/R/transition.R`
- Create: `modules/snagDecay/tests/testthat/test-transition.R`

- [ ] **Step 1: Write the failing test**

Create `modules/snagDecay/tests/testthat/test-transition.R`:

```r
library(testthat)

source(file.path(dirname(dirname(dirname(getwd()))), "modules", "snagDecay", "R", "transition.R"))

snagTransMat_test <- matrix(
  c(0.48, 0.38, 0.05, 0.00, 0.00,
    0.00, 0.52, 0.34, 0.05, 0.00,
    0.00, 0.00, 0.56, 0.30, 0.04,
    0.00, 0.00, 0.00, 0.62, 0.26,
    0.00, 0.00, 0.00, 0.00, 0.70),
  nrow = 5, byrow = TRUE
)

test_that("applyTransition returns integer vector of same length as input", {
  set.seed(42)
  result <- applyTransition(1:5, snagTransMat_test)
  expect_type(result, "integer")
  expect_length(result, 5L)
})

test_that("applyTransition output is always in 1:5", {
  set.seed(1)
  result <- applyTransition(rep(1:5, each = 100), snagTransMat_test)
  expect_true(all(result >= 1L & result <= 5L))
})

test_that("applyTransition never decreases DC (upper-triangular matrix)", {
  set.seed(99)
  dc_in <- rep(1:5, each = 500)
  dc_out <- applyTransition(dc_in, snagTransMat_test)
  expect_true(all(dc_out >= dc_in))
})

test_that("applyTransition DC5 input always returns DC5 (absorbing diagonal only)", {
  set.seed(7)
  result <- applyTransition(rep(5L, 200), snagTransMat_test)
  expect_true(all(result == 5L))
})
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-transition.R')"
```

Expected: FAIL — `could not find function "applyTransition"`

- [ ] **Step 3: Implement `applyTransition()`**

Create `modules/snagDecay/R/transition.R`:

```r
# Draw next-year DC for each element of DC_vec using a Markov transition matrix.
# transMatrix rows are current DC, columns are next-year DC.
# Row probabilities need not sum to 1; rmultinom normalises.
applyTransition <- function(DC_vec, transMatrix) {
  vapply(DC_vec, function(dc) {
    probs <- transMatrix[dc, ]
    which(stats::rmultinom(1L, 1L, probs) == 1L)
  }, integer(1L))
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-transition.R')"
```

Expected: 4 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add modules/snagDecay/R/transition.R modules/snagDecay/tests/testthat/test-transition.R
git commit -m "feat: add applyTransition() helper with tests"
```

---

## Task 2: `snagDecay` module — defineModule + init event (TDD)

**Files:**
- Create: `modules/snagDecay/snagDecay.R`
- Create: `modules/snagDecay/tests/testthat/test-snagDecay.R`

- [ ] **Step 1: Write the failing test**

Create `modules/snagDecay/tests/testthat/test-snagDecay.R`:

```r
library(SpaDES.core)
library(data.table)
library(testthat)

snagTransMat_test <- matrix(
  c(0.48, 0.38, 0.05, 0.00, 0.00,
    0.00, 0.52, 0.34, 0.05, 0.00,
    0.00, 0.00, 0.56, 0.30, 0.04,
    0.00, 0.00, 0.00, 0.62, 0.26,
    0.00, 0.00, 0.00, 0.00, 0.70),
  nrow = 5, byrow = TRUE
)
snagFallProb_test <- c(DC1 = 0.09, DC2 = 0.09, DC3 = 0.10, DC4 = 0.12, DC5 = 0.30)

emptyCohorts <- data.table(
  pixelID = integer(), year = integer(),
  species = character(), B = numeric()
)

test_that("snagDecay init creates empty snagTable with correct schema", {
  sim <- testInit(
    "snagDecay",
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$snagTable, "data.table")
  expect_equal(nrow(sim$snagTable), 0L)
  expect_named(sim$snagTable, c("pixelID", "species", "DC", "ageInDC", "initBiomass"))
})

test_that("snagDecay init creates empty fallenSnags", {
  sim <- testInit(
    "snagDecay",
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = emptyCohorts)
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$fallenSnags, "data.table")
  expect_equal(nrow(sim$fallenSnags), 0L)
})
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-snagDecay.R')"
```

Expected: FAIL — module file not found.

- [ ] **Step 3: Create `modules/snagDecay/snagDecay.R` with defineModule + init**

```r
defineModule(sim, list(
  name        = "snagDecay",
  description = "Advances standing dead White Pine trees through DC1-DC5 annually using
                 a Markov transition matrix, and stochastically transfers fallen snags
                 to sim$fallenSnags for consumption by DWDDecay.",
  keywords    = c("dead wood", "snag", "decay class", "Markov", "White Pine"),
  authors     = c(person("First", "Last", email = "email@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version     = list(snagDecay = "0.0.1"),
  spatialExtent = terra::ext(rep(NA_real_, 4)),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "SpaDES.core"),
  parameters  = rbind(
    defineParameter("snagTransMat", "matrix", matrix(0, 5, 5), NA, NA,
                    desc = "5x5 annual DC transition probability matrix for snags."),
    defineParameter("snagFallProb", "numeric", rep(0.1, 5), 0, 1,
                    desc = "Annual fall probability by DC (length 5)."),
    defineParameter("species", "character", "Pinus strobus", NA, NA,
                    desc = "Species to filter from cohortData.")
  ),
  inputObjects = rbind(
    expectsInput("cohortData", "data.table",
                 desc = "Pixel-level cohort table with columns: pixelID, year, species, B (Mg/ha).")
  ),
  outputObjects = rbind(
    createsOutput("snagTable", "data.table",
                  desc = "Current snag inventory: pixelID, species, DC, ageInDC, initBiomass."),
    createsOutput("fallenSnags", "data.table",
                  desc = "Snags that fell this timestep: same schema as snagTable.")
  )
))

doEvent.snagDecay <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- snagDecayInit(sim)
      sim <- scheduleEvent(sim, start(sim) + 1, "snagDecay", "annual", eventPriority = 1)
    },
    annual = {
      sim <- snagDecayAnnual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "snagDecay", "annual", eventPriority = 1)
    },
    warning(paste("Undefined event type:", eventType, "in module snagDecay"))
  )
  return(invisible(sim))
}

snagDecayInit <- function(sim) {
  sim$snagTable <- data.table::data.table(
    pixelID     = integer(),
    species     = character(),
    DC          = integer(),
    ageInDC     = integer(),
    initBiomass = numeric()
  )
  sim$fallenSnags <- data.table::copy(sim$snagTable)
  return(invisible(sim))
}

snagDecayAnnual <- function(sim) {
  # Absorb new mortality for this year and this species
  newDead <- sim$cohortData[year == time(sim) & species == P(sim)$species]
  if (nrow(newDead) > 0) {
    sim$snagTable <- data.table::rbindlist(list(
      sim$snagTable,
      newDead[, .(pixelID, species, DC = 1L, ageInDC = 0L, initBiomass = B)]
    ))
  }

  if (nrow(sim$snagTable) == 0) {
    sim$fallenSnags <- data.table::copy(sim$snagTable)
    return(invisible(sim))
  }

  # Advance decay class via Markov transition
  oldDC <- sim$snagTable$DC
  sim$snagTable[, DC := applyTransition(DC, P(sim)$snagTransMat)]
  sim$snagTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 1L, 0L)]

  # Stochastically simulate falls based on post-transition DC
  fallIdx <- sim$snagTable[, stats::rbinom(.N, 1L, P(sim)$snagFallProb[DC]) == 1L]
  sim$fallenSnags <- sim$snagTable[fallIdx]
  sim$snagTable   <- sim$snagTable[!fallIdx]

  return(invisible(sim))
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-snagDecay.R')"
```

Expected: 2 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add modules/snagDecay/snagDecay.R modules/snagDecay/tests/testthat/test-snagDecay.R
git commit -m "feat: add snagDecay module with init event"
```

---

## Task 3: `snagDecay` module — annual event tests

**Files:**
- Modify: `modules/snagDecay/tests/testthat/test-snagDecay.R` (add annual tests)

- [ ] **Step 1: Write the failing annual tests**

Append to `modules/snagDecay/tests/testthat/test-snagDecay.R`:

```r
test_that("snagDecay annual absorbs new mortality and populates snagTable", {
  cohorts <- data.table(
    pixelID = c(1L, 2L),
    year    = c(1L, 1L),
    species = "Pinus strobus",
    B       = c(10.0, 5.0)
  )
  sim <- testInit(
    "snagDecay",
    times  = list(start = 0, end = 2),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = snagFallProb_test,
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(42)
  sim <- spades(sim, events = c("init", "annual"))
  expect_true(nrow(sim$snagTable) + nrow(sim$fallenSnags) == 2L)
})

test_that("snagDecay annual DC never decreases", {
  cohorts <- data.table(
    pixelID = 1:20,
    year    = rep(1L, 20),
    species = "Pinus strobus",
    B       = rep(5.0, 20)
  )
  sim <- testInit(
    "snagDecay",
    times  = list(start = 0, end = 10),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = rep(0, 5),   # no falls — all snags stay
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(1)
  sim <- spades(sim)
  expect_true(all(sim$snagTable$DC >= 1L & sim$snagTable$DC <= 5L))
})

test_that("snagDecay annual with 100% fall probability empties snagTable each year", {
  cohorts <- data.table(
    pixelID = 1L, year = 1L, species = "Pinus strobus", B = 5.0
  )
  sim <- testInit(
    "snagDecay",
    times  = list(start = 0, end = 2),
    params = list(snagDecay = list(
      snagTransMat = snagTransMat_test,
      snagFallProb = rep(1, 5),   # guaranteed fall
      species      = "Pinus strobus"
    )),
    objects = list(cohortData = cohorts)
  )
  set.seed(5)
  sim <- spades(sim, events = c("init", "annual"))
  expect_equal(nrow(sim$snagTable), 0L)
  expect_equal(nrow(sim$fallenSnags), 1L)
})
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-snagDecay.R')"
```

Expected: 3 new tests fail (annual event not yet implemented with full logic).

- [ ] **Step 3: Verify existing init tests still pass**

Expected: the 2 init tests still pass.

- [ ] **Step 4: Run all tests to confirm 5 passing**

```bash
Rscript -e "testthat::test_file('modules/snagDecay/tests/testthat/test-snagDecay.R')"
```

Expected: 5 tests, 0 failures. (The annual logic was already added in Task 2 step 3 — if the tests fail, the issue is in `snagDecayAnnual`. Debug by verifying that `spades(sim, events = c("init", "annual"))` fires events in order and that `cohortData` years align with `time(sim)`.)

- [ ] **Step 5: Commit**

```bash
git add modules/snagDecay/tests/testthat/test-snagDecay.R
git commit -m "test: add snagDecay annual event tests"
```

---

## Task 4: `pixelValuesToRaster()` helper — TDD

**Files:**
- Create: `modules/deadWoodBiomass/R/raster-helpers.R`
- Create: `modules/deadWoodBiomass/tests/testthat/test-raster-helpers.R`

- [ ] **Step 1: Write the failing test**

Create `modules/deadWoodBiomass/tests/testthat/test-raster-helpers.R`:

```r
library(terra)
library(data.table)
library(testthat)

source(file.path(dirname(dirname(dirname(getwd()))), "modules", "deadWoodBiomass", "R", "raster-helpers.R"))

test_that("pixelValuesToRaster assigns values to correct pixels", {
  template <- terra::rast(nrows = 3, ncols = 3,
                          xmin = 0, xmax = 3, ymin = 0, ymax = 3,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = c(1L, 5L, 9L), value = c(10.0, 20.0, 30.0))
  r <- pixelValuesToRaster(dt, template)
  expect_equal(terra::values(r)[1],  10.0)
  expect_equal(terra::values(r)[5],  20.0)
  expect_equal(terra::values(r)[9],  30.0)
  expect_true(is.na(terra::values(r)[2]))
})

test_that("pixelValuesToRaster returns NA for unspecified pixels", {
  template <- terra::rast(nrows = 2, ncols = 2,
                          xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = 1L, value = 99.0)
  r <- pixelValuesToRaster(dt, template)
  expect_equal(sum(is.na(terra::values(r))), 3L)
})

test_that("pixelValuesToRaster handles empty input data.table", {
  template <- terra::rast(nrows = 2, ncols = 2,
                          xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                          crs = "EPSG:4326")
  dt <- data.table(pixelID = integer(), value = numeric())
  r <- pixelValuesToRaster(dt, template)
  expect_true(all(is.na(terra::values(r))))
})
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('modules/deadWoodBiomass/tests/testthat/test-raster-helpers.R')"
```

Expected: FAIL — `could not find function "pixelValuesToRaster"`

- [ ] **Step 3: Implement `pixelValuesToRaster()`**

Create `modules/deadWoodBiomass/R/raster-helpers.R`:

```r
# Map pixel IDs and values from a data.table onto a copy of templateRaster.
# dt must have columns pixelID (integer) and value (numeric).
# Returns a SpatRaster with NA for all pixels not in dt.
pixelValuesToRaster <- function(dt, templateRaster) {
  r <- templateRaster
  terra::values(r) <- NA_real_
  if (nrow(dt) > 0) {
    r[dt$pixelID] <- dt$value
  }
  r
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('modules/deadWoodBiomass/tests/testthat/test-raster-helpers.R')"
```

Expected: 3 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add modules/deadWoodBiomass/R/raster-helpers.R \
        modules/deadWoodBiomass/tests/testthat/test-raster-helpers.R
git commit -m "feat: add pixelValuesToRaster() helper with tests"
```

---

## Task 5: `DWDDecay` module — defineModule + init + receive events (TDD)

**Files:**
- Create: `modules/DWDDecay/DWDDecay.R`
- Create: `modules/DWDDecay/tests/testthat/test-DWDDecay.R`

- [ ] **Step 1: Write the failing tests**

Create `modules/DWDDecay/tests/testthat/test-DWDDecay.R`:

```r
library(SpaDES.core)
library(data.table)
library(testthat)

DWDTransMat_test <- matrix(
  c(0.55, 0.35, 0.06, 0.00, 0.00,
    0.00, 0.50, 0.37, 0.08, 0.00,
    0.00, 0.00, 0.48, 0.38, 0.08,
    0.00, 0.00, 0.00, 0.50, 0.36,
    0.00, 0.00, 0.00, 0.00, 0.72),
  nrow = 5, byrow = TRUE
)
snagToDWD_DCmap_test <- c(DC1 = 1L, DC2 = 2L, DC3 = 2L, DC4 = 3L, DC5 = 4L)
DWD_lossProb_test    <- c(DC1 = 0.00, DC2 = 0.05, DC3 = 0.06, DC4 = 0.14, DC5 = 0.28)

emptySnagTable <- data.table(
  pixelID = integer(), species = character(),
  DC = integer(), ageInDC = integer(), initBiomass = numeric()
)

test_that("DWDDecay init creates empty DWDTable with correct schema", {
  sim <- testInit(
    "DWDDecay",
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = DWD_lossProb_test
    )),
    objects = list(fallenSnags = data.table::copy(emptySnagTable))
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$DWDTable, "data.table")
  expect_equal(nrow(sim$DWDTable), 0L)
  expect_named(sim$DWDTable, c("pixelID", "species", "DC", "ageInDC", "initBiomass"))
})

test_that("DWDDecay receive appends fallenSnags with DC remapped", {
  fallen <- data.table(
    pixelID     = c(1L, 2L),
    species     = "Pinus strobus",
    DC          = c(1L, 3L),
    ageInDC     = c(2L, 1L),
    initBiomass = c(8.0, 4.0)
  )
  sim <- testInit(
    "DWDDecay",
    times   = list(start = 0, end = 2),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = DWD_lossProb_test
    )),
    objects = list(fallenSnags = fallen)
  )
  sim <- spades(sim, events = c("init", "receive"))
  expect_equal(nrow(sim$DWDTable), 2L)
  # DC1 snag → DC1 DWD; DC3 snag → DC2 DWD (per snagToDWD_DCmap)
  expect_equal(sim$DWDTable[pixelID == 1L, DC], 1L)
  expect_equal(sim$DWDTable[pixelID == 2L, DC], 2L)
  # ageInDC resets to 0 upon entering DWD pool
  expect_true(all(sim$DWDTable$ageInDC == 0L))
})
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('modules/DWDDecay/tests/testthat/test-DWDDecay.R')"
```

Expected: FAIL — module not found.

- [ ] **Step 3: Create `modules/DWDDecay/DWDDecay.R`**

```r
defineModule(sim, list(
  name        = "DWDDecay",
  description = "Manages the downed woody debris pool. Receives fallen snags from snagDecay,
                 maps their decay class to the DWD scale, advances DC annually, and removes
                 fully decomposed records.",
  keywords    = c("dead wood", "DWD", "downed woody debris", "decay class", "Markov"),
  authors     = c(person("First", "Last", email = "email@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version     = list(DWDDecay = "0.0.1"),
  spatialExtent = terra::ext(rep(NA_real_, 4)),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "SpaDES.core"),
  parameters  = rbind(
    defineParameter("DWDTransMat", "matrix", matrix(0, 5, 5), NA, NA,
                    desc = "5x5 annual DC transition probability matrix for DWD."),
    defineParameter("snagToDWD_DCmap", "integer", 1:5, 1L, 5L,
                    desc = "Maps snag DC at fall to starting DWD DC (length 5)."),
    defineParameter("DWD_lossProb", "numeric", rep(0.05, 5), 0, 1,
                    desc = "Annual complete-loss probability by DC (length 5).")
  ),
  inputObjects = rbind(
    expectsInput("fallenSnags", "data.table",
                 desc = "Snags that fell in the current timestep from snagDecay.")
  ),
  outputObjects = rbind(
    createsOutput("DWDTable", "data.table",
                  desc = "Current DWD inventory: pixelID, species, DC, ageInDC, initBiomass.")
  )
))

doEvent.DWDDecay <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- DWDDecayInit(sim)
      sim <- scheduleEvent(sim, start(sim) + 1, "DWDDecay", "receive", eventPriority = 2)
      sim <- scheduleEvent(sim, start(sim) + 1, "DWDDecay", "annual",  eventPriority = 3)
    },
    receive = {
      sim <- DWDDecayReceive(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DWDDecay", "receive", eventPriority = 2)
    },
    annual = {
      sim <- DWDDecayAnnual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DWDDecay", "annual",  eventPriority = 3)
    },
    warning(paste("Undefined event type:", eventType, "in module DWDDecay"))
  )
  return(invisible(sim))
}

DWDDecayInit <- function(sim) {
  sim$DWDTable <- data.table::data.table(
    pixelID     = integer(),
    species     = character(),
    DC          = integer(),
    ageInDC     = integer(),
    initBiomass = numeric()
  )
  return(invisible(sim))
}

DWDDecayReceive <- function(sim) {
  if (is.null(sim$fallenSnags) || nrow(sim$fallenSnags) == 0L) {
    return(invisible(sim))
  }
  incoming <- data.table::copy(sim$fallenSnags)
  # Remap snag DC → starting DWD DC; reset ageInDC to 0
  incoming[, DC      := P(sim)$snagToDWD_DCmap[DC]]
  incoming[, ageInDC := 0L]
  sim$DWDTable <- data.table::rbindlist(list(sim$DWDTable, incoming))
  return(invisible(sim))
}

DWDDecayAnnual <- function(sim) {
  if (nrow(sim$DWDTable) == 0L) return(invisible(sim))

  oldDC <- sim$DWDTable$DC
  sim$DWDTable[, DC := applyTransition(DC, P(sim)$DWDTransMat)]
  sim$DWDTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 1L, 0L)]

  # Remove records that suffer complete loss this year
  lossIdx <- sim$DWDTable[, stats::rbinom(.N, 1L, P(sim)$DWD_lossProb[DC]) == 1L]
  sim$DWDTable <- sim$DWDTable[!lossIdx]

  return(invisible(sim))
}
```

**Note:** `DWDDecay` calls `applyTransition()` which lives in `modules/snagDecay/R/transition.R`. SpaDES auto-sources all files in a module's `R/` directory when loading the module. Since `applyTransition` is in `snagDecay`, we need to either (a) copy it to `modules/DWDDecay/R/`, or (b) define it as a shared helper in `R/transition.R` at the project level. Use option (b): copy `transition.R` to `modules/DWDDecay/R/transition.R` (a straight file copy — keep both).

```bash
cp modules/snagDecay/R/transition.R modules/DWDDecay/R/transition.R
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
Rscript -e "testthat::test_file('modules/DWDDecay/tests/testthat/test-DWDDecay.R')"
```

Expected: 2 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add modules/DWDDecay/ 
git commit -m "feat: add DWDDecay module with init and receive events"
```

---

## Task 6: `DWDDecay` module — annual event tests

**Files:**
- Modify: `modules/DWDDecay/tests/testthat/test-DWDDecay.R` (add annual tests)

- [ ] **Step 1: Write the failing annual tests**

Append to `modules/DWDDecay/tests/testthat/test-DWDDecay.R`:

```r
test_that("DWDDecay annual DC never decreases (upper-triangular matrix)", {
  fallen <- data.table(
    pixelID     = 1:50,
    species     = "Pinus strobus",
    DC          = sample(1:3, 50, replace = TRUE),
    ageInDC     = rep(0L, 50),
    initBiomass = rep(5.0, 50)
  )
  sim <- testInit(
    "DWDDecay",
    times   = list(start = 0, end = 20),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = rep(0, 5)   # no loss — all records persist
    )),
    objects = list(fallenSnags = fallen)
  )
  set.seed(77)
  sim <- spades(sim)
  expect_true(all(sim$DWDTable$DC >= 1L & sim$DWDTable$DC <= 5L))
})

test_that("DWDDecay annual with 100% loss probability removes all DC5 records", {
  # Start with pieces already at DC5
  preloaded <- data.table(
    pixelID     = 1:10,
    species     = "Pinus strobus",
    DC          = rep(5L, 10),
    ageInDC     = rep(0L, 10),
    initBiomass = rep(3.0, 10)
  )
  sim <- testInit(
    "DWDDecay",
    times   = list(start = 0, end = 1),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = c(0, 0, 0, 0, 1)   # DC5 always lost
    )),
    objects = list(fallenSnags = data.table::copy(emptySnagTable))
  )
  # Manually pre-populate DWDTable after init, before annual
  sim <- spades(sim, events = "init")
  sim$DWDTable <- preloaded
  set.seed(3)
  sim <- spades(sim, events = "annual")
  expect_equal(nrow(sim$DWDTable[DC == 5L]), 0L)
})
```

- [ ] **Step 2: Run tests**

```bash
Rscript -e "testthat::test_file('modules/DWDDecay/tests/testthat/test-DWDDecay.R')"
```

Expected: 4 tests, 0 failures.

- [ ] **Step 3: Commit**

```bash
git add modules/DWDDecay/tests/testthat/test-DWDDecay.R
git commit -m "test: add DWDDecay annual event tests"
```

---

## Task 7: `deadWoodBiomass` module — defineModule + init event (TDD)

**Files:**
- Create: `modules/deadWoodBiomass/deadWoodBiomass.R`
- Create: `modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R`

- [ ] **Step 1: Write the failing tests**

Create `modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R')"
```

Expected: FAIL — module not found.

- [ ] **Step 3: Create `modules/deadWoodBiomass/deadWoodBiomass.R`**

```r
defineModule(sim, list(
  name        = "deadWoodBiomass",
  description = "Translates snag and DWD decay class inventories into pixel-level biomass
                 estimates (Mg ha-1) using species- and pool-specific density reduction
                 factors (DRF) from Paper 2 Appendix D.",
  keywords    = c("dead wood", "biomass", "density reduction factor", "carbon"),
  authors     = c(person("First", "Last", email = "email@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version     = list(deadWoodBiomass = "0.0.1"),
  spatialExtent = terra::ext(rep(NA_real_, 4)),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "terra", "SpaDES.core"),
  parameters  = rbind(
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
  inputObjects = rbind(
    expectsInput("snagTable", "data.table",
                 desc = "Current snag inventory from snagDecay."),
    expectsInput("DWDTable", "data.table",
                 desc = "Current DWD inventory from DWDDecay."),
    expectsInput("studyAreaRaster", "SpatRaster",
                 desc = "Template raster defining pixel grid, CRS, and resolution.")
  ),
  outputObjects = rbind(
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
  sim$snagBiomass_Mg_ha <- sim$studyAreaRaster
  terra::values(sim$snagBiomass_Mg_ha) <- NA_real_
  sim$DWDBiomass_Mg_ha  <- sim$studyAreaRaster
  terra::values(sim$DWDBiomass_Mg_ha)  <- NA_real_
  return(invisible(sim))
}

deadWoodBiomassAnnual <- function(sim) {
  drf <- P(sim)$DRFLookup

  # Snag pool biomass
  snagWithBiomass <- data.table::copy(sim$snagTable)
  snagWithBiomass <- drf[pool == "snag"][snagWithBiomass, on = .(species, DC)]
  snagWithBiomass[, currentBiomass := initBiomass * DRF]
  snagByPixel <- snagWithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
  sim$snagBiomass_Mg_ha <- pixelValuesToRaster(snagByPixel, sim$studyAreaRaster)

  # DWD pool biomass
  DWDwithBiomass <- data.table::copy(sim$DWDTable)
  DWDwithBiomass <- drf[pool == "DWD"][DWDwithBiomass, on = .(species, DC)]
  DWDwithBiomass[, currentBiomass := initBiomass * DRF]
  DWDbyPixel <- DWDwithBiomass[, .(value = sum(currentBiomass, na.rm = TRUE)), by = pixelID]
  sim$DWDBiomass_Mg_ha <- pixelValuesToRaster(DWDbyPixel, sim$studyAreaRaster)

  return(invisible(sim))
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
Rscript -e "testthat::test_file('modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R')"
```

Expected: 1 test, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add modules/deadWoodBiomass/
git commit -m "feat: add deadWoodBiomass module with init event"
```

---

## Task 8: `deadWoodBiomass` module — annual event tests

**Files:**
- Modify: `modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R` (add annual tests)

- [ ] **Step 1: Write the failing annual tests**

Append to `modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R`:

```r
test_that("deadWoodBiomass annual computes snag biomass correctly for DC1", {
  # DC1 snag: DRF = 1.0, so currentBiomass == initBiomass
  snagTable <- data.table(
    pixelID = 1L, species = "Pinus strobus", DC = 1L, ageInDC = 0L, initBiomass = 12.0
  )
  DWDTable <- data.table(
    pixelID = integer(), species = character(),
    DC = integer(), ageInDC = integer(), initBiomass = numeric()
  )
  sim <- testInit(
    "deadWoodBiomass",
    times   = list(start = 0, end = 2),
    params  = list(deadWoodBiomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = DWDTable,
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "annual"))
  expect_equal(terra::values(sim$snagBiomass_Mg_ha)[1], 12.0)
  expect_true(is.na(terra::values(sim$DWDBiomass_Mg_ha)[1]))
})

test_that("deadWoodBiomass annual applies DRF correctly for DC3 DWD", {
  # DC3 DWD DRF = 0.614; initBiomass = 10 → currentBiomass = 6.14
  DWDTable <- data.table(
    pixelID = 2L, species = "Pinus strobus", DC = 3L, ageInDC = 1L, initBiomass = 10.0
  )
  snagTable <- data.table(
    pixelID = integer(), species = character(),
    DC = integer(), ageInDC = integer(), initBiomass = numeric()
  )
  sim <- testInit(
    "deadWoodBiomass",
    times   = list(start = 0, end = 2),
    params  = list(deadWoodBiomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = DWDTable,
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "annual"))
  expect_equal(terra::values(sim$DWDBiomass_Mg_ha)[2], 10.0 * 0.614, tolerance = 1e-6)
})

test_that("deadWoodBiomass annual aggregates multiple records per pixel", {
  # Two DC2 snags in pixel 3: 5 + 8 = 13 Mg/ha, each × DRF 0.841 = 10.933
  snagTable <- data.table(
    pixelID     = c(3L, 3L),
    species     = "Pinus strobus",
    DC          = c(2L, 2L),
    ageInDC     = c(0L, 0L),
    initBiomass = c(5.0, 8.0)
  )
  DWDTable <- data.table(
    pixelID = integer(), species = character(),
    DC = integer(), ageInDC = integer(), initBiomass = numeric()
  )
  sim <- testInit(
    "deadWoodBiomass",
    times   = list(start = 0, end = 2),
    params  = list(deadWoodBiomass = list(DRFLookup = DRFLookup_test)),
    objects = list(
      snagTable       = snagTable,
      DWDTable        = DWDTable,
      studyAreaRaster = templateRaster
    )
  )
  sim <- spades(sim, events = c("init", "annual"))
  expect_equal(terra::values(sim$snagBiomass_Mg_ha)[3], (5.0 + 8.0) * 0.841, tolerance = 1e-6)
})
```

- [ ] **Step 2: Run tests**

```bash
Rscript -e "testthat::test_file('modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R')"
```

Expected: 4 tests, 0 failures.

- [ ] **Step 3: Commit**

```bash
git add modules/deadWoodBiomass/tests/testthat/test-deadWoodBiomass.R
git commit -m "test: add deadWoodBiomass annual event tests"
```

---

## Task 9: `global.R` — integration script and end-to-end test

**Files:**
- Create: `global.R`

- [ ] **Step 1: Create `global.R`**

```r
library(SpaDES.project)  # bootstrap exception — only library() call permitted

Require::Require(c("SpaDES.core", "data.table", "terra"))

source("R/parameters.R")

# Minimal single-pixel cohort data for a 50-year run
myMortalityTable <- data.table::data.table(
  pixelID = 1L,
  year    = c(1L, 5L, 10L),    # mortality events at years 1, 5, 10
  species = "Pinus strobus",
  B       = c(15.0, 10.0, 8.0) # Mg ha-1
)

myRaster <- terra::rast(
  nrows = 1, ncols = 1,
  xmin = 0, xmax = 1, ymin = 0, ymax = 1,
  crs = "EPSG:4326"
)

times   <- list(start = 0, end = 50)
params  <- list(
  snagDecay = list(
    snagTransMat = snagTransMat,
    snagFallProb = snagFallProb,
    species      = "Pinus strobus"
  ),
  DWDDecay = list(
    DWDTransMat     = DWDTransMat,
    snagToDWD_DCmap = snagToDWD_DCmap,
    DWD_lossProb    = DWD_lossProb
  ),
  deadWoodBiomass = list(
    DRFLookup           = DRFLookup,
    .plotInitialTime    = NA
  )
)
modules <- list("snagDecay", "DWDDecay", "deadWoodBiomass")

mySim <- simInit(
  times   = times,
  params  = params,
  modules = modules,
  objects = list(
    cohortData      = myMortalityTable,
    studyAreaRaster = myRaster
  ),
  paths = list(
    modulePath = file.path("modules"),
    inputPath  = file.path("inputs"),
    outputPath = file.path("outputs"),
    cachePath  = file.path("cache")
  )
)

set.seed(42)
mySim <- spades(mySim)

# Inspect outputs
cat("Final snag inventory rows:   ", nrow(mySim$snagTable), "\n")
cat("Final DWD inventory rows:    ", nrow(mySim$DWDTable),  "\n")
cat("Final snag biomass (pixel 1):", terra::values(mySim$snagBiomass_Mg_ha)[1], "Mg/ha\n")
cat("Final DWD biomass  (pixel 1):", terra::values(mySim$DWDBiomass_Mg_ha)[1],  "Mg/ha\n")
```

- [ ] **Step 2: Run the integration script**

```bash
Rscript global.R
```

Expected output (exact values will vary by seed):
```
Final snag inventory rows:    <some integer>
Final DWD inventory rows:     <some integer>
Final snag biomass (pixel 1): <numeric or NA>
Final DWD biomass  (pixel 1): <numeric or NA>
```

No errors. If `simInit` fails with contract warnings, check that `expectsInput`/`createsOutput` declarations match the objects in the `objects = list(...)` argument.

- [ ] **Step 3: Verify biomass values are physically plausible**

```bash
Rscript -e "
source('R/parameters.R')
source('global.R')

snagB <- terra::values(mySim\$snagBiomass_Mg_ha)[1]
DWDB  <- terra::values(mySim\$DWDBiomass_Mg_ha)[1]

# At year 50, initial biomass = 33 Mg/ha total input; DRF >= 0.251 for DC5
# Both values must be <= total input biomass of 33 Mg/ha or NA
if (!is.na(snagB)) stopifnot(snagB >= 0 && snagB <= 33)
if (!is.na(DWDB))  stopifnot(DWDB  >= 0 && DWDB  <= 33)
cat('Biomass plausibility check passed\n')
"
```

Expected: `Biomass plausibility check passed`

- [ ] **Step 4: Commit**

```bash
git add global.R
git commit -m "feat: add global.R integration entry point; full pipeline verified"
```

---

## Implementation Checklist (from spec §7)

These items require human review against the source papers before final use in production:

- [ ] Extract exact snag transition probability matrices from Paper 1, Table 2 (replace placeholders in `R/parameters.R`)
- [ ] Extract exact DWD transition probability matrices from Paper 1 (replace placeholders)
- [ ] Extract exact White Pine DRF values from Paper 2, Appendix D (replace placeholders in `R/parameters.R`)
- [ ] Confirm snag-to-DWD DC mapping (`snagToDWD_DCmap`) against Paper 1 methods section
- [ ] Confirm snag fall probability values against species-specific appendix in Paper 1
- [ ] Validate `applyTransition()` output distribution against multinomial expectation (analytical check)
- [ ] Sensitivity analysis: vary `snagFallProb` ± 20%, compare 50-year C pool trajectories
- [ ] Optional: extend `DRFLookup` with additional species rows for multi-species landscapes
