# CLAUDE.md

## Project Overview

This is a **Stata/Mata demonstration repository** that showcases advanced fixed-effects regression techniques using the `reghdfe` package. The single script serves as a reference implementation for manually performing fixed-effects regression via Mata (Stata's matrix programming language), including iterative solvers and optimization methods.

## Repository Structure

```
atheodor/
├── mata-demo-update.do   # Main Stata/Mata demonstration script
└── CLAUDE.md             # This file
```

This is a single-file project with no subdirectories, build system, or configuration files.

## Technology Stack

- **Language:** Stata `.do` script with embedded Mata code blocks
- **Dependencies:**
  - Stata (econometric software)
  - `reghdfe` package (high-dimensional fixed effects regression)
  - `reghdfe.mata` library (Mata implementation of FixedEffects class)

## What the Code Does

`mata-demo-update.do` demonstrates a complete workflow:

1. **Setup** - Loads the built-in `auto` dataset and defines regression variables
   - Dependent variable: `price`
   - Independent variables: `weight`, `gear`
   - Fixed-effect (absorb) variables: `turn`, `trunk`
2. **Standard regression** - Runs `reghdfe` for comparison
3. **Manual Mata implementation** - Steps through the same regression manually:
   - Creates a `FixedEffects` object and configures absorb variables
   - Partials out fixed effects using `partial_out()`
   - Loads data into Mata matrices via `st_data()`
   - Applies symmetric Kaczmarz transformation with conjugate gradient acceleration
   - Standardizes data before solving
   - Solves OLS via QR decomposition (`qrsolve`)
   - Unstandardizes coefficients using stored standard deviations

## Stata/Mata Conventions Used

- **Local macros** use backtick-quote syntax: `` `depvar' ``, `` `indepvars' ``
- **Mata blocks** are delimited by `mata:` ... `end`
- **Object-oriented Mata** uses dot notation: `HDFE.absvars`, `HDFE.solution.stdevs`
- **Matrix slicing**: `data[., 1]` (column 1), `data[., 2::k]` (columns 2 through k)
- **Function pointers** assigned with `&`: `fun_transform = &transform_sym_kaczmarz()`
- **Comments**: `//` style in Mata blocks

## Development Workflow

There is no build, test, or CI/CD pipeline. To run the script:

1. Open Stata
2. Execute: `do mata-demo-update.do`

The script requires the `reghdfe` package to be installed in Stata (`ssc install reghdfe`).

## Notes for AI Assistants

- This is an **educational/reference** project, not a production software system
- There are no tests, linters, formatters, or build steps to run
- The `reghdfe.mata` file is loaded at runtime via `include` with `adopath` — it is not included in this repository
- When modifying Stata code, preserve the local macro syntax (backtick-quote pairs must be balanced)
- The order of operations in the Mata block matters — `partial_out()` must run before manual data manipulation to ensure `HDFE.solution` is properly initialized
- Standardization via `standardize_data()` must be undone after solving to get correct coefficient estimates
