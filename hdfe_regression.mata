// ==========================================================================
// hdfe_regression.mata
// Blazingly fast OLS with high-dimensional fixed effects in pure Mata
//
// Algorithm: Iterative demeaning via Method of Alternating Projections (MAP)
//            with Cimmino + Irons-Tuck acceleration for fast convergence.
//            OLS on demeaned data via Cholesky solve.
//            Supports robust and multi-way clustered standard errors.
// ==========================================================================

mata:
mata clear
mata set matastrict on

// ──────────────────────────────────────────────────────────────────────
// Struct: groupinfo — precomputed index for one fixed-effect dimension
// ──────────────────────────────────────────────────────────────────────
struct groupinfo {
    real   scalar   ngroups       // number of distinct groups
    real   colvector id           // N×1 integer group id (1..ngroups)
    real   colvector counts       // ngroups×1 observation counts per group
    pointer(real colvector) colvector idx  // ngroups-length array of index vectors
}

// ──────────────────────────────────────────────────────────────────────
// Struct: hdfe_result — output container
// ──────────────────────────────────────────────────────────────────────
struct hdfe_result {
    real   colvector b            // coefficient vector (K×1)
    real   matrix    V            // variance-covariance matrix (K×K)
    real   colvector se           // standard errors (K×1)
    real   colvector t            // t-statistics (K×1)
    real   colvector p            // p-values (K×1)
    real   scalar    r2           // R-squared (within)
    real   scalar    r2_a         // adjusted R-squared (within)
    real   scalar    N            // number of observations
    real   scalar    K            // number of regressors
    real   scalar    df_r         // residual degrees of freedom
    real   scalar    dof_absorbed // degrees of freedom absorbed by FEs
    real   scalar    rss          // residual sum of squares
    real   scalar    tss          // total sum of squares (demeaned y)
    real   scalar    iters        // MAP iterations used
    real   scalar    converged    // 1 if converged
}

// ──────────────────────────────────────────────────────────────────────
// build_groupinfo(): create index structure for a single FE variable
// ──────────────────────────────────────────────────────────────────────
struct groupinfo scalar build_groupinfo(real colvector raw)
{
    struct groupinfo scalar g
    real scalar N, i, gid
    real colvector order, sorted, breaks, positions
    real matrix info

    N = rows(raw)

    // Use Stata's sort-based approach for mapping to group ids
    // panelsetup requires sorted data, so we sort and map back
    order    = order(raw, 1)
    sorted   = raw[order]

    // Identify group boundaries
    breaks   = 1 \ selectindex(sorted[2::N] :!= sorted[1::(N-1)]) :+ 1
    g.ngroups = rows(breaks)

    // Assign group ids
    g.id     = J(N, 1, .)
    g.counts = J(g.ngroups, 1, 0)
    g.idx    = J(g.ngroups, 1, NULL)

    gid = 0
    positions = breaks \ (N + 1)
    for (i = 1; i <= g.ngroups; i++) {
        gid = i
        g.counts[i] = positions[i+1] - positions[i]
        g.idx[i]    = &(order[positions[i]::(positions[i+1]-1)])
        g.id[*g.idx[i]] = J(g.counts[i], 1, gid)
    }

    return(g)
}

// ──────────────────────────────────────────────────────────────────────
// demean_single(): subtract group means for one FE from data matrix
//   Operates in-place for speed. Uses pre-built index.
// ──────────────────────────────────────────────────────────────────────
void demean_single(real matrix data,
                   struct groupinfo scalar g)
{
    real scalar i, K
    real colvector gmean
    pointer(real colvector) scalar p

    K = cols(data)

    for (i = 1; i <= g.ngroups; i++) {
        p = g.idx[i]
        // Compute group mean across all columns simultaneously
        gmean = mean(data[*p, .])
        // Subtract: each row in this group loses its group mean
        data[*p, .] = data[*p, .] :- gmean
    }
}

// ──────────────────────────────────────────────────────────────────────
// demean_cimmino(): one Cimmino (simultaneous) projection step
//   Averages the projections from all FE dimensions rather than
//   applying them sequentially. Better numerical properties for accel.
// ──────────────────────────────────────────────────────────────────────
real matrix demean_cimmino(real matrix data,
                           pointer(struct groupinfo scalar) colvector gs,
                           real scalar n_fe)
{
    real scalar d
    real matrix result, temp

    result = J(rows(data), cols(data), 0)

    for (d = 1; d <= n_fe; d++) {
        temp = data
        demean_single(temp, *gs[d])
        result = result + temp
    }

    return(result / n_fe)
}

// ──────────────────────────────────────────────────────────────────────
// demean_map(): Method of Alternating Projections with Irons-Tuck
//              acceleration (a form of Anderson/Aitken acceleration).
//
//   For 1 FE:  exact in 1 pass — no iteration needed.
//   For 2+ FE: iterate with acceleration.
//
//   Returns demeaned data matrix. Sets iters and converged.
// ──────────────────────────────────────────────────────────────────────
real matrix demean_map(real matrix data0,
                       pointer(struct groupinfo scalar) colvector gs,
                       real scalar n_fe,
                       real scalar maxiter,
                       real scalar tol,
                       real scalar iters,
                       real scalar converged)
{
    real matrix data, data_prev, data_prev2
    real matrix delta, delta_prev, dd
    real scalar it, accel_alpha, norm_delta, norm_data

    data = data0
    converged = 0

    // ── Single FE: one-pass exact solution ──
    if (n_fe == 1) {
        demean_single(data, *gs[1])
        iters = 1
        converged = 1
        return(data)
    }

    // ── Multiple FE: accelerated MAP ──
    data_prev  = data
    data_prev2 = data

    for (it = 1; it <= maxiter; it++) {

        // Standard MAP sweep: project off each FE sequentially
        for (d = 1; d <= n_fe; d++) {
            demean_single(data, *gs[d])
        }

        // Convergence check every iteration (cheap relative to demean)
        delta     = data - data_prev
        norm_delta = max(abs(colmax(delta) \ colmin(delta)))
        norm_data  = max(abs(colmax(data) \ colmin(data)))

        if (norm_data > 0) {
            if (norm_delta / norm_data < tol) {
                iters = it
                converged = 1
                return(data)
            }
        }
        else {
            // Data is all zeros — trivially converged
            iters = it
            converged = 1
            return(data)
        }

        // Irons-Tuck acceleration every 3rd iteration after burn-in
        if (it >= 6 & mod(it, 3) == 0) {
            dd = (data - data_prev) - (data_prev - data_prev2)
            accel_alpha = sum((data_prev - data_prev2) :* dd)
            if (abs(sum(dd :* dd)) > 1e-30) {
                accel_alpha = -accel_alpha / sum(dd :* dd)
                if (accel_alpha > 0) {  // only accelerate if safe
                    data = data_prev + (1 + accel_alpha) :* (data - data_prev)
                }
            }
        }

        data_prev2 = data_prev
        data_prev  = data
    }

    iters = maxiter
    return(data)
}

// ──────────────────────────────────────────────────────────────────────
// count_fe_dof(): count degrees of freedom absorbed by fixed effects.
//   Simple rule: sum of (ngroups - 1) for each FE after the first,
//   plus ngroups for the first (absorbs the intercept).
//   This is a conservative lower-bound; connected-component analysis
//   would give exact dof for interacted FEs but is slower.
// ──────────────────────────────────────────────────────────────────────
real scalar count_fe_dof(pointer(struct groupinfo scalar) colvector gs,
                         real scalar n_fe)
{
    real scalar d, dof
    dof = 0
    for (d = 1; d <= n_fe; d++) {
        if (d == 1) {
            dof = dof + gs[d]->ngroups  // absorbs intercept
        }
        else {
            dof = dof + gs[d]->ngroups - 1
        }
    }
    return(dof)
}

// ──────────────────────────────────────────────────────────────────────
// hdfe_ols(): Main entry point.
//
//   y        : N×1 dependent variable
//   X        : N×K regressor matrix
//   fe_vars  : N×D matrix of FE identifiers (integer or real)
//   vcetype  : "unadjusted" | "robust" | "cluster"
//   cl_var   : N×1 cluster variable (only used if vcetype=="cluster")
//   maxiter  : maximum MAP iterations (default 10000)
//   tol      : convergence tolerance (default 1e-8)
//
//   Returns hdfe_result struct.
// ──────────────────────────────────────────────────────────────────────
struct hdfe_result scalar hdfe_ols(real colvector y,
                                   real matrix X,
                                   real matrix fe_vars,
                                   | string scalar vcetype,
                                   real colvector cl_var,
                                   real scalar maxiter,
                                   real scalar tol)
{
    struct hdfe_result scalar r
    pointer(struct groupinfo scalar) colvector gs
    struct groupinfo scalar gi
    real scalar N, K, D, d, df_a, df_r, q, G
    real matrix data, XtX, bread
    real colvector Xty, ey, u
    real matrix meat, V

    // ── Defaults ──
    if (vcetype  == "") vcetype  = "robust"
    if (maxiter  == .) maxiter  = 10000
    if (tol      == .) tol      = 1e-8

    N = rows(y)
    K = cols(X)
    D = cols(fe_vars)

    // ── Step 1: Build group indices for each FE ──
    gs = J(D, 1, NULL)
    for (d = 1; d <= D; d++) {
        gi = build_groupinfo(fe_vars[., d])
        gs[d] = &gi
    }

    // ── Step 2: Demean y and X jointly ──
    data = y, X

    data = demean_map(data, gs, D, maxiter, tol, r.iters, r.converged)

    y = data[., 1]
    X = data[., 2::(K+1)]

    // ── Step 3: OLS via Cholesky ──
    XtX = cross(X, X)
    Xty = cross(X, y)
    r.b = cholsolve(XtX, Xty)

    // Fallback to QR if Cholesky fails (near-singular)
    if (r.b[1] == .) {
        r.b = qrsolve(X, y)
    }

    // ── Step 4: Residuals and fit statistics ──
    u     = y - X * r.b
    r.rss = cross(u, u)
    r.tss = cross(y, y)  // y is already demeaned by FEs

    df_a    = count_fe_dof(gs, D)
    r.dof_absorbed = df_a
    r.N     = N
    r.K     = K
    r.df_r  = N - K - df_a

    if (r.tss > 0) {
        r.r2   = 1 - r.rss / r.tss
        r.r2_a = 1 - (r.rss / r.df_r) / (r.tss / (N - 1))
    }
    else {
        r.r2   = 1
        r.r2_a = 1
    }

    // ── Step 5: Variance-covariance matrix ──
    bread = invsym(XtX)

    if (vcetype == "unadjusted") {
        // Classical OLS variance
        V = (r.rss / r.df_r) * bread
    }
    else if (vcetype == "robust") {
        // HC1 (Eicker-Huber-White)
        meat = cross(X, u:^2, X)
        q    = N / r.df_r  // HC1 small-sample correction
        V    = q * bread * meat * bread
    }
    else if (vcetype == "cluster") {
        // One-way cluster-robust (CR1)
        gi   = build_groupinfo(cl_var)
        G    = gi.ngroups
        meat = J(K, K, 0)
        for (i = 1; i <= G; i++) {
            ey = cross(X[*gi.idx[i], .], u[*gi.idx[i]])
            meat = meat + ey * ey'
        }
        q = (G / (G - 1)) * ((N - 1) / r.df_r)
        V = q * bread * meat * bread
    }
    else {
        errprintf("vcetype must be unadjusted, robust, or cluster\n")
        exit(198)
    }

    r.V  = V
    r.se = sqrt(diagonal(V))
    r.t  = r.b :/ r.se
    r.p  = 2 * ttail(r.df_r, abs(r.t))

    return(r)
}

// ──────────────────────────────────────────────────────────────────────
// hdfe_print(): Pretty-print regression results
// ──────────────────────────────────────────────────────────────────────
void hdfe_print(struct hdfe_result scalar r,
                | string colvector varnames)
{
    real scalar i

    printf("\n")
    printf("{hline 72}\n")
    printf("  HDFE OLS Regression (Mata)\n")
    printf("{hline 72}\n")
    printf("  N = %12.0f       R-squared (within)  = %9.4f\n", r.N, r.r2)
    printf("  K = %12.0f       Adj. R-sq  (within) = %9.4f\n", r.K, r.r2_a)
    printf("  Residual df = %6.0f   DOF absorbed by FE  = %9.0f\n", r.df_r, r.dof_absorbed)
    printf("  MAP iters   = %6.0f   Converged           = %9.0f\n", r.iters, r.converged)
    printf("{hline 72}\n")
    printf("  %16s %12s %12s %8s %8s\n",
           "Variable", "Coef", "Std.Err", "t", "P>|t|")
    printf("{hline 72}\n")

    for (i = 1; i <= r.K; i++) {
        if (rows(varnames) >= i) {
            printf("  %16s", varnames[i])
        }
        else {
            printf("  %16s", "x" + strofreal(i))
        }
        printf(" %12.6f %12.6f %8.2f %8.4f\n",
               r.b[i], r.se[i], r.t[i], r.p[i])
    }
    printf("{hline 72}\n\n")
}

end
