// ==========================================================================
// hdfe_demo.do
// Demo: blazingly fast high-dimensional FE regression in pure Mata
// ==========================================================================
clear all
cls
set more off

// ──────────────────────────────────────────────────────────────────────
// 1) Load the Mata library
// ──────────────────────────────────────────────────────────────────────
do hdfe_regression.mata

// ══════════════════════════════════════════════════════════════════════
// EXAMPLE 1: Small data — auto dataset, two-way FE
// ══════════════════════════════════════════════════════════════════════
sysuse auto, clear

mata:
    // Pull data into Mata
    y  = st_data(., "price")
    X  = st_data(., "weight gear_ratio")
    fe = st_data(., "turn trunk")

    // Run HDFE OLS with robust SE (default)
    res = hdfe_ols(y, X, fe)
    hdfe_print(res, ("weight" \ "gear_ratio"))

    // Now with classical (unadjusted) SE
    res2 = hdfe_ols(y, X, fe, "unadjusted")
    hdfe_print(res2, ("weight" \ "gear_ratio"))
end

// Compare with reghdfe for validation
reghdfe price weight gear_ratio, a(turn trunk) vce(robust)

// ══════════════════════════════════════════════════════════════════════
// EXAMPLE 2: Large simulated data with 3 high-dimensional FEs
// ══════════════════════════════════════════════════════════════════════
clear
set seed 42
local N = 500000

set obs `N'

// Generate 3 high-dimensional FEs
gen long fe1 = ceil(runiform() * 5000)
gen long fe2 = ceil(runiform() * 3000)
gen long fe3 = ceil(runiform() * 1000)

// True DGP: y = 2*x1 - 0.5*x2 + fe1_effect + fe2_effect + fe3_effect + noise
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double y  = 2 * x1 - 0.5 * x2 + 0.1 * fe1 + 0.05 * fe2 + 0.01 * fe3 + rnormal(0, 2)

// Time the Mata implementation
timer clear 1
timer on 1

mata:
    y  = st_data(., "y")
    X  = st_data(., "x1 x2")
    fe = st_data(., "fe1 fe2 fe3")

    res = hdfe_ols(y, X, fe, "robust", ., 50000, 1e-8)
    hdfe_print(res, ("x1" \ "x2"))
end

timer off 1
display "Mata HDFE time:"
timer list 1

// Time reghdfe for comparison
timer clear 2
timer on 2
reghdfe y x1 x2, a(fe1 fe2 fe3) vce(robust)
timer off 2
display "reghdfe time:"
timer list 2

// ══════════════════════════════════════════════════════════════════════
// EXAMPLE 3: Clustered standard errors
// ══════════════════════════════════════════════════════════════════════
mata:
    cl = st_data(., "fe1")  // cluster on fe1
    res_cl = hdfe_ols(y, X, fe, "cluster", cl)
    hdfe_print(res_cl, ("x1" \ "x2"))
end

reghdfe y x1 x2, a(fe1 fe2 fe3) vce(cluster fe1)

exit
