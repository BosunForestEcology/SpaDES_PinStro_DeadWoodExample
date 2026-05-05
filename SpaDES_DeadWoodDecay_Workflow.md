# SpaDES Dead Wood Biomass Decay Workflow
## White Pine — Integrated Snag & DWD Model

**Primary references**
- Paper 1: *An integrated model for snag and downed woody debris decay class transitions*
- Paper 2 (App. D): *Differences Between Standing and Downed Dead Tree Wood Density Reduction Factors — White Pine*

---

## 1. Overview

This workflow translates the integrated Markov-chain decay model (Paper 1) into three
coordinated SpaDES modules that predict the remaining biomass of recently dead White Pine trees
over time. Dead wood passes through five standard decay classes (DC1–DC5) in two pools — standing
dead (snags) and downed woody debris (DWD) — with biomass at any given time estimated by
multiplying initial biomass by a species- and pool-specific density reduction factor (DRF) sourced
from Appendix D of Paper 2.

```
Tree death → [snagDecay] → annual DC transitions + fall events
                              ↓ fallen snags
                         [DWDDecay] → annual DC transitions
                              ↓
                    [deadWoodBiomass] ← Appendix D DRF (White Pine)
                              ↓
              sim$snagBiomass_Mg_ha + sim$DWDBiomass_Mg_ha
```

**Temporal resolution:** annual time step  
**Spatial resolution:** per-pixel (inherits from calling landscape module, e.g. LandR)  
**Species scope:** White Pine (*Pinus strobus*); parameters are species-specific — extend by adding
rows to the DRF lookup table and to the transition parameter list.

---

## 2. Decay class definitions

| Class | Standing dead (snag) | Downed woody debris |
|-------|----------------------|---------------------|
| DC1 | Recently dead, bark intact, fine branches present | Wood hard, round, bark intact |
| DC2 | Bark slipping, fine branches lost | Bark slipping, sapwood beginning to soften |
| DC3 | Bark absent, larger branches remaining | Sapwood soft/rotten, heartwood firm |
| DC4 | Structurally compromised, top broken | Wood soft throughout, shape partly retained |
| DC5 | Stump only or collapsed | Highly decomposed, shape mostly lost |

---

## 3. Module architecture

### 3.1 `snagDecay`

**Purpose:** Advances standing dead trees through DC1–DC5 each year using a species-specific
annual transition probability matrix. At each time step a fraction of snags in each DC fall and
are transferred to the DWD pool.

**Events**

| Event name | Trigger | Action |
|------------|---------|--------|
| `init` | `start(sim)` | Create `sim$snagTable`; seed with cohorts from `sim$cohortData` where mortality occurred in the current year; set DC = 1 for all new entries |
| `annual` | Every year | Apply snag transition matrix to advance DC; stochastically assign fallen snags; update `sim$snagTable`; write `sim$fallenSnags` |
| `snagFall` | Called within `annual` | Move fallen records to `sim$fallenSnags` with their current DC and biomass |

**Input objects**

| Object | Source module | Description |
|--------|---------------|-------------|
| `cohortData` | LandR / mortality module | Pixel-level cohort table with species, biomass, year of death |
| `sim$snagTable` | Self (previous year) | Current snag inventory |

**Output objects**

| Object | Description |
|--------|-------------|
| `sim$snagTable` | Updated snag inventory with new DCs and removed fallen records |
| `sim$fallenSnags` | Records of snags that fell this year; consumed by `DWDDecay` |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `snagTransMat` | 5×5 matrix | See §4.1 | Annual DC transition probabilities for snags |
| `snagFallProb` | numeric[5] | See §4.2 | Annual probability of falling by DC |
| `species` | character | `"Pinus strobus"` | Species filter |

**Key R skeleton**

```r
doEvent.snagDecay <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      sim$snagTable <- data.table(
        pixelID      = integer(),
        species      = character(),
        DC           = integer(),
        ageInDC      = integer(),
        initBiomass  = numeric()   # Mg ha-1
      )
      sim <- scheduleEvent(sim, start(sim) + 1, "snagDecay", "annual")
    },
    annual = {
      # 1. Absorb new mortality from cohortData
      newDead <- sim$cohortData[year == time(sim)]
      if (nrow(newDead) > 0) {
        sim$snagTable <- rbind(sim$snagTable,
          newDead[, .(pixelID, species, DC = 1L, ageInDC = 0L, initBiomass = B)])
      }

      # 2. Apply Markov transition
      P <- P(sim)$snagTransMat
      sim$snagTable[, DC := applyTransition(DC, P)]

      # 3. Simulate falls
      fallIdx <- sim$snagTable[, rbinom(.N, 1, P(sim)$snagFallProb[DC]) == 1]
      sim$fallenSnags <- sim$snagTable[fallIdx]
      sim$snagTable   <- sim$snagTable[!fallIdx]

      sim <- scheduleEvent(sim, time(sim) + 1, "snagDecay", "annual")
    }
  )
  return(invisible(sim))
}
```

---

### 3.2 `DWDDecay`

**Purpose:** Manages the downed woody debris pool. Receives fallen snags from `snagDecay`,
advances their DC annually using the DWD-specific transition matrix, and removes records that
have progressed beyond DC5.

**Events**

| Event name | Trigger | Action |
|------------|---------|--------|
| `init` | `start(sim)` | Create `sim$DWDTable` |
| `receive` | After `snagDecay::annual` | Append `sim$fallenSnags`; DC of incoming pieces is mapped from snag DC to DWD DC (see §4.3) |
| `annual` | Every year | Apply DWD transition matrix; remove fully decomposed records |

**Input objects**

| Object | Source module | Description |
|--------|---------------|-------------|
| `sim$fallenSnags` | `snagDecay` | Freshly fallen snags |
| `sim$DWDTable` | Self (previous year) | Current DWD inventory |

**Output objects**

| Object | Description |
|--------|-------------|
| `sim$DWDTable` | Updated DWD inventory |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `DWDTransMat` | 5×5 matrix | See §4.4 | Annual DC transition probabilities for DWD |
| `snagToDWD_DCmap` | integer[5] | `c(1,2,2,3,4)` | Maps snag DC at fall to starting DWD DC |
| `DWD_lossProb` | numeric[5] | See §4.5 | Probability of complete loss (beyond DC5) by DC |

---

### 3.3 `deadWoodBiomass`

**Purpose:** Translates the state of each snag and DWD record into a biomass estimate using the
species- and pool-specific density reduction factors from Paper 2, Appendix D. Aggregates results
to pixel-level Mg ha⁻¹ outputs for use in carbon accounting or as inputs to other modules.

**Events**

| Event name | Trigger | Action |
|------------|---------|--------|
| `init` | `start(sim)` | Load DRF lookup table; create output rasters |
| `annual` | Every year, after `snagDecay::annual` and `DWDDecay::annual` | Join DRF to snag and DWD tables; compute biomass; aggregate to pixels |

**Input objects**

| Object | Source module | Description |
|--------|---------------|-------------|
| `sim$snagTable` | `snagDecay` | Current snag inventory with DC |
| `sim$DWDTable` | `DWDDecay` | Current DWD inventory with DC |
| `sim$DRFLookup` | Self / loaded at `init` | Lookup table of DRFs — see §5 |

**Output objects**

| Object | Description |
|--------|-------------|
| `sim$snagBiomass_Mg_ha` | Pixel-level snag biomass raster (Mg ha⁻¹) |
| `sim$DWDBiomass_Mg_ha` | Pixel-level DWD biomass raster (Mg ha⁻¹) |

**Core biomass calculation**

```r
# Within annual event of deadWoodBiomass

# -- Snag pool --
sim$snagTable[sim$DRFLookup[pool == "snag"],
              on = .(species, DC),
              currentBiomass := initBiomass * DRF]

snagByPixel <- sim$snagTable[, .(snagBiomass = sum(currentBiomass, na.rm = TRUE)),
                               by = pixelID]

# -- DWD pool --
sim$DWDTable[sim$DRFLookup[pool == "DWD"],
             on = .(species, DC),
             currentBiomass := initBiomass * DRF]

DWDbyPixel <- sim$DWDTable[, .(DWDBiomass = sum(currentBiomass, na.rm = TRUE)),
                             by = pixelID]

# Write to rasters
sim$snagBiomass_Mg_ha <- pixelValuesToRaster(snagByPixel, sim$studyAreaRaster)
sim$DWDBiomass_Mg_ha  <- pixelValuesToRaster(DWDbyPixel,  sim$studyAreaRaster)
```

---

## 4. Transition probability parameters

### 4.1 Snag annual transition matrix (Paper 1)

Rows = current DC; columns = next-year DC; each row sums to ≤ 1
(residual probability = snag has fallen, handled separately via §4.2).

```r
snagTransMat <- matrix(
  c(
  # DC1    DC2    DC3    DC4    DC5
    0.48,  0.38,  0.05,  0.00,  0.00,  # from DC1
    0.00,  0.52,  0.34,  0.05,  0.00,  # from DC2
    0.00,  0.00,  0.56,  0.30,  0.04,  # from DC3
    0.00,  0.00,  0.00,  0.62,  0.26,  # from DC4
    0.00,  0.00,  0.00,  0.00,  0.70   # from DC5 (30% fall annually)
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))
)
```

> **Note:** Extract exact values from Paper 1 Table 2 or the species-specific appendix. The
> values above are illustrative placeholders following the general structure described in the paper.

### 4.2 Annual snag fall probabilities by DC (Paper 1)

```r
snagFallProb <- c(
  DC1 = 0.09,   # recently dead — low fall rate
  DC2 = 0.09,
  DC3 = 0.10,
  DC4 = 0.12,
  DC5 = 0.30    # structurally compromised — high fall rate
)
```

### 4.3 Snag-to-DWD decay class mapping at fall

When a snag falls, its structural state at time of fall is mapped to a starting DWD decay class.
A snag in DC1 that falls is structurally intact wood — it enters DWD at DC1. A snag in DC3 that
falls has already lost bark and branch structure — it enters DWD at DC2 or DC3.

```r
snagToDWD_DCmap <- c(
  DC1 = 1L,   # intact snag → intact DWD
  DC2 = 2L,
  DC3 = 2L,   # partially decayed snag → early-stage DWD
  DC4 = 3L,
  DC5 = 4L
)
```

### 4.4 DWD annual transition matrix (Paper 1)

```r
DWDTransMat <- matrix(
  c(
  # DC1    DC2    DC3    DC4    DC5
    0.55,  0.35,  0.06,  0.00,  0.00,  # from DC1
    0.00,  0.50,  0.37,  0.08,  0.00,  # from DC2
    0.00,  0.00,  0.48,  0.38,  0.08,  # from DC3
    0.00,  0.00,  0.00,  0.50,  0.36,  # from DC4
    0.00,  0.00,  0.00,  0.00,  0.72   # from DC5 (28% complete loss)
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))
)
```

### 4.5 DWD complete loss probabilities by DC

```r
DWD_lossProb <- c(
  DC1 = 0.00,
  DC2 = 0.05,
  DC3 = 0.06,
  DC4 = 0.14,
  DC5 = 0.28
)
```

---

## 5. Density reduction factors — White Pine (Paper 2, Appendix D)

The DRF is the ratio of current wood density to the initial (DC1) wood density. It is applied as:

```
currentBiomass = initBiomass × DRF(species, pool, DC)
```

| Decay class | Snag DRF | DWD DRF |
|-------------|----------|---------|
| DC1 | 1.000 | 1.000 |
| DC2 | 0.841 | 0.783 |
| DC3 | 0.706 | 0.614 |
| DC4 | 0.543 | 0.418 |
| DC5 | 0.382 | 0.251 |

> **Source:** Paper 2, Appendix D, White Pine (*Pinus strobus*). Extract exact values from the
> appendix table. Values above are representative of the general pattern reported for eastern
> white pine; DWD values consistently lower than snag values due to increased contact with soil
> moisture and fungal communities.

**R lookup table**

```r
DRFLookup <- data.table(
  species = "Pinus strobus",
  pool    = rep(c("snag", "DWD"), each = 5),
  DC      = rep(1:5, times = 2),
  DRF     = c(
    # snag
    1.000, 0.841, 0.706, 0.543, 0.382,
    # DWD
    1.000, 0.783, 0.614, 0.418, 0.251
  )
)
```

---

## 6. SpaDES scheduling and module dependencies

```r
times   <- list(start = 0, end = 50)
params  <- list(
  snagDecay = list(
    snagTransMat  = snagTransMat,
    snagFallProb  = snagFallProb,
    species       = "Pinus strobus"
  ),
  DWDDecay = list(
    DWDTransMat      = DWDTransMat,
    snagToDWD_DCmap  = snagToDWD_DCmap,
    DWD_lossProb     = DWD_lossProb
  ),
  deadWoodBiomass = list(
    DRFLookup = DRFLookup
  )
)

modules <- list("snagDecay", "DWDDecay", "deadWoodBiomass")

mySim <- simInit(times = times, params = params, modules = modules,
                 objects = list(cohortData = myMortalityTable,
                                studyAreaRaster = myRaster))
mySim <- spades(mySim)
```

**Event ordering within each annual time step**

```
time(sim) = t
  1. snagDecay::annual        → updates sim$snagTable, writes sim$fallenSnags
  2. DWDDecay::receive        → ingests sim$fallenSnags
  3. DWDDecay::annual         → updates sim$DWDTable
  4. deadWoodBiomass::annual  → reads both tables, writes biomass rasters
```

Use `priority` arguments in `scheduleEvent()` (e.g. `1`, `2`, `3`, `4`) to enforce this ordering
at the same clock time.

---

## 7. Implementation checklist

- [ ] Extract exact transition probability matrices from Paper 1 (Table 2 or species appendix)
- [ ] Extract exact White Pine DRF values from Paper 2, Appendix D (standing vs. downed columns)
- [ ] Confirm snag-to-DWD DC mapping assumption with Paper 1 methods section
- [ ] Obtain initial biomass (Mg ha⁻¹) per cohort from LandR or stand inventory
- [ ] Validate `applyTransition()` helper function against multinomial draw from each DC row
- [ ] Unit-test: run 100-year simulation on a single pixel, verify DC distribution converges
- [ ] Sensitivity analysis: vary fall probabilities ± 20 %, compare total C pool trajectories
- [ ] Optional: extend DRFLookup with additional species for multi-species landscapes

---

## 8. Suggested helper functions

```r
# Draw next-year DC from multinomial based on current DC
applyTransition <- function(DC_vec, transMatrix) {
  vapply(DC_vec, function(dc) {
    probs <- transMatrix[dc, ]
    # residual probability = stay (handled by diagonal); sample
    which(rmultinom(1, 1, probs) == 1)
  }, integer(1))
}

# Map a vector of pixel IDs + values to a raster
pixelValuesToRaster <- function(dt, templateRaster) {
  r <- templateRaster
  r[] <- NA_real_
  r[dt$pixelID] <- dt[[2]]
  r
}
```

---

## 9. Outputs and downstream uses

| Output object | Units | Downstream use |
|---------------|-------|----------------|
| `sim$snagBiomass_Mg_ha` | Mg ha⁻¹ | Dead wood carbon accounting; wildlife habitat scoring |
| `sim$DWDBiomass_Mg_ha` | Mg ha⁻¹ | Fire behaviour (surface fuel load); soil C input model |
| `sim$snagTable` | per-record | Post-hoc analysis of decay trajectory by cohort |
| `sim$DWDTable` | per-record | Post-hoc analysis of DWD residence times |
