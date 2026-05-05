# Module Context Template

Use this template to create a `CLAUDE.md` file inside a specific SpaDES module directory.
Copy and fill in the sections below. Delete sections that don't apply.

---

## Module Identity

**Name:** `moduleName`
**Purpose:** One sentence — what ecological process does this module simulate and why does it exist in this pipeline?
**Position in pipeline:** What modules run before this one? What modules consume its outputs?

---

## Reference Docs

Before editing this module, read:
- `.claude/spades.ai/core/module-anatomy.md` — `defineModule()` structure, `doEvent()` patterns
- `.claude/spades.ai/core/simlist-events.md` — `simList` accessors, event scheduling
- `.claude/spades.ai/core/caching-reproducibility.md` — `Cache()`, `prepInputs()`, `suppliedElsewhere()`

---

## Inputs

| Object | Class | Provided by | Description |
|--------|-------|-------------|-------------|
| `studyArea` | `SpatVector` | user / upstream module | Study area polygon |
| `rasterToMatch` | `SpatRaster` | upstream module | Template raster for CRS/resolution |
| `objectName` | `data.table` | `upstreamModule` | What this object contains |

---

## Outputs

| Object | Class | Consumed by | Description |
|--------|-------|-------------|-------------|
| `outputName` | `data.table` | `downstreamModule` | What this object contains |

---

## Events

| Event | When it fires | What it does |
|-------|---------------|-------------|
| `init` | Once at `start(sim)` | Downloads inputs, creates initial state |
| `decay` | Every year | Applies decay equations to cohort data |
| `plot` | Per `.plotInterval` | Visualises current state |

---

## Non-obvious Implementation Details

Document gotchas here — things an AI or new developer would likely get wrong:

- **Key invariant:** `sim$decayCohorts` is keyed on `(pixelGroup, speciesCode)` — any merge must preserve this key or downstream modules silently produce wrong results.
- **Why X is done this way:** Brief explanation of a non-obvious design choice.
- **Temporary workaround:** Description of workaround + link to tracking issue.

---

## Known Dependencies / Load Order

List any required load order constraints:

- Must initialise after `Biomass_core` (needs `sim$cohortData` to exist)
- Must initialise before `Biomass_validation` (provides `sim$decayOutput`)

---

## Active Issues

| Issue | Status | Workaround |
|-------|--------|------------|
| Link to issue tracker | open/in-progress | Brief description of workaround |

---

## Resource Links

- Module repo: `https://github.com/org/moduleName`
- Issue tracker: `https://github.com/org/moduleName/issues`
- Related modules: list names and repos
