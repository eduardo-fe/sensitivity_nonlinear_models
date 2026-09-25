*===============================================================================
*  SENSITIVITY ANALYSIS FOR NONLINEAR MODELS — Version 18
*
*  Based on V17 validate-then-bisect design. Changes from V17:
*
*  1. Genuine adjacency: crossings require consecutively successful
*     grid points. A failed evaluation breaks the bracket.
*
*  2. Sign-aware crossing detection for explain-away: when two adjacent
*     scan points are both significant but the ME sign has flipped,
*     the |t| function must have crossed z_alpha within that interval.
*
*  3. Explicit anchor check at rho_D = 0: if the criterion is already
*     met at zero confounding, the threshold is reported as 0.
*
*  4. Failure handling: failed bisection midpoints shrink toward the
*     known-good (criterion-not-met) side.
*
*  5. Multinomial benchmarking: sens_abc uses the main-model r_o,
*     not the reduced-model r_bench.
*
*  PRESERVED from V17: sign-change first, explain-away capped at
*  sign-change threshold via rdbound(). damage_dir logic unchanged.
*
*===============================================================================

clear all
set more off
version 16


capture program drop sens_abc
program define sens_abc
    syntax , rd(real) ry(real) rho(real)

    global sens_feasible = 0
    global sens_a        = .
    global sens_b        = .
    global sens_c        = .

    local denom = 1 - (`rho')^2
    if abs(`denom') < 1e-8 exit

    local av  = (`rd' - `ry'*`rho') / `denom'
    local bv  = (`ry' - `rd'*`rho') / `denom'
    local csq = 1 - `av'^2 - `bv'^2 - 2*`av'*`bv'*`rho'
    if `csq' < -1e-10 exit
    local csq = max(0, `csq')

    global sens_a        = `av'
    global sens_b        = `bv'
    global sens_c        = sqrt(`csq')
    global sens_feasible = 1
end


*===============================================================================
* sens_build_z
*   Constructs _Z_hyp_ = a*Dresid + b*Yresid + c*eps.
*   Uses globals: sens_a  sens_b  sens_c
*===============================================================================
capture program drop sens_build_z
program define sens_build_z
    syntax , dresid(varname) yresid(varname) touse(varname) [ epsvar(varname) ]

    capture drop _Z_hyp_

    * Use fixed epsilon: explicit option > global $sens_eps > fresh draw
    if "`epsvar'" != "" {
        local eps_to_use "`epsvar'"
    }
    else if "$sens_eps" != "" {
        local eps_to_use "$sens_eps"
    }
    else {
        tempvar eps_fresh
        qui gen double `eps_fresh' = rnormal() if `touse'
        local eps_to_use "`eps_fresh'"
    }

    qui gen double _Z_hyp_ = $sens_a * `dresid'      ///
                            + $sens_b * `yresid'      ///
                            + $sens_c * `eps_to_use'  if `touse'
end


*===============================================================================
* sens_fit
*   Fits model with _Z_hyp_ and retrieves AME via margins.
*   Uses global: sens_X
*   Globals set: sens_me  sens_se  sens_failed
*===============================================================================
capture program drop sens_fit
program define sens_fit
    syntax , model(string) yvar(varname) dvar(varname) touse(varname) ///
             atm(integer) baseoutcome(integer) [ outcome(string) ]

    global sens_me     = .
    global sens_se     = .
    global sens_failed = 0

    capture {
        if      "`model'" == "logit"   qui logit   `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "probit"  qui probit  `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "cloglog" qui cloglog `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "poisson" qui poisson `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "nbreg"   qui nbreg   `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "ologit"  qui ologit  `yvar' `dvar' $sens_X _Z_hyp_ if `touse'
        else if "`model'" == "mlogit" {
            * iterate(50) bounds the wall-clock cost of any single fit. Near
            * the feasibility boundary, _Z_hyp_ becomes nearly collinear with
            * Dr/Yr, which can put mlogit's likelihood on a "nearly flat
            * section that stretches out to infinity" (Stata's own FAQ on ML
            * convergence) coefficients march toward huge values while the
            * log-likelihood barely moves, so the optimizer can run for a very
            * long time without ever satisfying the convergence criterion.
            * Capping iterations turns "occasionally takes minutes" into
            * "fails fast", which sens_failed below already knows how to
            * handle like any other infeasible point. e(converged) must be
            * checked explicitly because hitting the iterate() cap still
            * leaves _rc==0 -- Stata returns a results table either way.
            qui mlogit `yvar' `dvar' $sens_X _Z_hyp_ if `touse', ///
                baseoutcome(`baseoutcome') iterate(50)
            if e(converged) == 0 error 430
        }

        if `atm' == 1 {
            if "`outcome'" != "" {
                qui margins, dydx(`dvar') atmeans predict(outcome(`outcome'))
            }
            else qui margins, dydx(`dvar') atmeans
        }
        else {
            if "`outcome'" != "" {
                qui margins, dydx(`dvar') predict(outcome(`outcome'))
            }
            else qui margins, dydx(`dvar')
        }

        global sens_me    = r(table)[1,1]
        global sens_se    = r(table)[2,1]
        capture global sens_delta = _b[_Z_hyp_]
    }

    if _rc != 0 | $sens_me == . | $sens_se == . global sens_failed = 1
    if $sens_failed == 1 global sens_delta = .
end


*===============================================================================
* sens_bisect   
*
*   Drop-in replacement for the original sens_bisect. Same call
*   signature, same globals set on exit (sens_rdstar, sens_mestar,
*   sens_tstar, sens_bstatus), PLUS one new global:
*     sens_bmatch -- 1 if a single crossing was found and it fell on
*                    the side the "damaging direction" heuristic
*                    predicted; 0 if it fell on the other side;
*                    missing if zero or multiple crossings were found.
*
*   Instead of jumping straight to the feasibility boundary and
*   bisecting based on the analytic damage_dir prediction, this version:
*     1. Computes the feasible rho_D range exactly as before (capped to
*        [0, rdbound] on the damage_dir-predicted side when rdbound()
*        is supplied, exactly mirroring the original's behaviour).
*     2. Scans nscan points across that range, recording success/failure
*        and whether the criterion is met at each successful point.
*     3. Locates where the criterion's truth value changes between
*        ADJACENT successful points -- the empirical crossing(s).
*     4. If exactly one crossing is found, bisects within that narrow,
*        already-verified bracket (both endpoints known to succeed),
*        which also means bisection is never asked to step into a
*        region the scan hasn't already shown to be numerically stable.
*     5. If zero crossings: status "robust". If more than one: status
*        "non-monotone: N crossings", reporting the threshold closest
*        to rho_D=0 (the most conservative / smallest-magnitude one).
*===============================================================================
capture program drop sens_bisect
program define sens_bisect
    syntax , ry(real) rho(real) meshort(real) model(string)      ///
             yvar(varname) dvar(varname) touse(varname)           ///
             dresid(varname) yresid(varname) zval(real)           ///
             criterion(string) atm(integer) baseoutcome(integer)  ///
             [ outcome(string) rdbound(real 99) nscan(integer 21) ///
               scanmargin(real 0.001) VERBose ]

    local bisect_tol = 0.0001
    local sign0       = sign(`meshort')

    global sens_rdstar  = .
    global sens_mestar  = .
    global sens_tstar   = .
    global sens_bstatus = "unknown"
    global sens_bmatch  = .

    *---------------------------------------------------------------------------
    * Feasible bracket for rho_D, given fixed ry and rho (same algebra
    * as the original).
    *---------------------------------------------------------------------------
    local disc_feas = (1 - (`rho')^2) * (1 - (`ry')^2)
    if `disc_feas' < 0 {
        global sens_bstatus = "no feasible region"
        exit
    }

    local rd_lo_feas = max(-0.999, `rho'*`ry' - sqrt(`disc_feas'))
    local rd_hi_feas = min( 0.999, `rho'*`ry' + sqrt(`disc_feas'))

    if `rd_hi_feas' <= `rd_lo_feas' {
        global sens_bstatus = "no feasible region"
        exit
    }

    sens_abc, rd(`rd_hi_feas') ry(`ry') rho(`rho')
    if $sens_feasible == 0 {
        local rd_hi_feas = `rd_hi_feas' - 0.001
    }
    sens_abc, rd(`rd_lo_feas') ry(`ry') rho(`rho')
    if $sens_feasible == 0 {
        local rd_lo_feas = `rd_lo_feas' + 0.001
    }
    if `rd_hi_feas' <= `rd_lo_feas' {
        global sens_bstatus = "no feasible region"
        exit
    }

    *---------------------------------------------------------------------------
    * damage_dir: computed both as a diagnostic label AND, when rdbound()
    * is supplied, to know which side to cap -- exactly mirroring the
    * original's use of rdbound() (cap the damaging side at rdbound,
    * treat rho_D=0 as the safe anchor on the other side).
    *---------------------------------------------------------------------------
    local damage_dir = sign(`sign0' * `ry')
    if `damage_dir' == 0 local damage_dir = 1

    if `rdbound' != 99 {
        if `damage_dir' > 0 {
            local rd_hi_feas = min(`rd_hi_feas', `rdbound')
            local rd_lo_feas = 0
        }
        else {
            local rd_lo_feas = max(`rd_lo_feas', `rdbound')
            local rd_hi_feas = 0
        }
        if `rd_hi_feas' <= `rd_lo_feas' {
            global sens_bstatus = "no feasible region"
            exit
        }
    }

    *---------------------------------------------------------------------------
    * Pull the scan range in from the TRUE feasibility edges by scanmargin
    * (a fraction of the feasible width), to avoid the degenerate region
    * where c^2 -> 0 and _Z_hyp_ becomes nearly a deterministic linear
    * combination of Dr and Yr -- exactly the regime where mlogit's
    * optimizer can stall or fail to converge cleanly. This is NOT applied
    * to a rho_D=0 anchor side (0 is always safe, no confounding at all),
    * only to whichever end(s) are still sitting at the true feasibility
    * boundary after the rdbound() capping above.
    *---------------------------------------------------------------------------
    local feas_width = `rd_hi_feas' - `rd_lo_feas'
    if `rd_lo_feas' != 0 {
        local rd_lo_feas = `rd_lo_feas' + `scanmargin' * `feas_width'
    }
    if `rd_hi_feas' != 0 {
        local rd_hi_feas = `rd_hi_feas' - `scanmargin' * `feas_width'
    }
    if `rd_hi_feas' <= `rd_lo_feas' {
        global sens_bstatus = "no feasible region"
        exit
    }

    *---------------------------------------------------------------------------
    * Coarse scan across the feasible range determined above.
    *---------------------------------------------------------------------------
    tempname scan
    matrix `scan' = J(`nscan', 4, .)
    * col 1 = rho_D, col 2 = success (0/1), col 3 = criterion met (0/1),
    * col 4 = sign of ME (V18: used for sign-aware crossing detection)

    forval s = 0/`=`nscan'-1' {
        local rd_s = `rd_lo_feas' + `s' * (`rd_hi_feas' - `rd_lo_feas') / (`nscan' - 1)

        if "`verbose'" != "" {
            di as text "  [sens_bisect scan] model=`model' outcome=`outcome' " ///
                "criterion=`criterion' ry=" %5.2f `ry' " scan_pt=" `s'+1 "/`nscan'" ///
                " rho_D=" %6.3f `rd_s'
        }

        sens_abc, rd(`rd_s') ry(`ry') rho(`rho')
        if $sens_feasible == 0 {
            matrix `scan'[`s'+1, 1] = `rd_s'
            matrix `scan'[`s'+1, 2] = 0
            continue
        }

        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')

        timer clear 99
        timer on 99
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        timer off 99
        qui timer list 99
        if "`verbose'" != "" {
            di as text "      -> fit took " %5.2f r(t99) " sec"
        }

        capture drop _Z_hyp_

        matrix `scan'[`s'+1, 1] = `rd_s'

        if $sens_failed == 0 & $sens_me != . {
            matrix `scan'[`s'+1, 2] = 1
            local tstat_s = $sens_me / $sens_se
            local crit_s  = 0
            if "`criterion'" == "insignif" & abs(`tstat_s') <  `zval'  local crit_s = 1
            if "`criterion'" == "signflip" & sign($sens_me) != `sign0' local crit_s = 1
            matrix `scan'[`s'+1, 3] = `crit_s'
            matrix `scan'[`s'+1, 4] = sign($sens_me)
        }
        else {
            matrix `scan'[`s'+1, 2] = 0
        }
    }

    *---------------------------------------------------------------------------
    * V18 FIX 3: Explicit anchor check at rho_D = 0.
    * If the criterion is already met at zero confounding, the threshold
    * IS zero. Evaluated before crossing detection.
    *---------------------------------------------------------------------------
    local anchor_crit = 0
    sens_abc, rd(0) ry(`ry') rho(`rho')
    if $sens_feasible == 1 {
        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        capture drop _Z_hyp_
        if $sens_failed == 0 & $sens_me != . {
            local tstat_0 = $sens_me / $sens_se
            if "`criterion'" == "insignif" & abs(`tstat_0') <  `zval'  local anchor_crit = 1
            if "`criterion'" == "signflip" & sign($sens_me) != `sign0' local anchor_crit = 1
        }
    }

    *---------------------------------------------------------------------------
    * Find crossings among successfully-evaluated, ordered points.
    * V18 FIX 1: genuine adjacency — a failed grid point breaks the bracket.
    * V18 FIX 2: sign-aware crossing for explain-away — when two adjacent
    *   points are both significant (crit=0) but the ME sign has flipped,
    *   the |t| function must have crossed z_alpha within that interval.
    *---------------------------------------------------------------------------
    local n_success  = 0
    local n_crossings = 0
    local first_cross_lo = .
    local first_cross_hi = .
    local cross_closest_to_zero = .
    local cross_lo_at_closest   = .
    local cross_hi_at_closest   = .

    local prev_rd   = .
    local prev_crit = .
    local prev_sign = .
    local prev_ok   = 0

    forval s = 1/`nscan' {
        local ok_s = `scan'[`s', 2]
        if `ok_s' == 1 {
            local n_success = `n_success' + 1
            local rd_s   = `scan'[`s', 1]
            local crit_s = `scan'[`s', 3]
            local sign_s = `scan'[`s', 4]

            * V18 FIX 1: only form bracket if PREVIOUS grid point was
            * also successful (genuinely adjacent)
            if `prev_ok' == 1 {
                local is_crossing = 0

                * Standard crossing: criterion status changes
                if `crit_s' != `prev_crit' {
                    local is_crossing = 1
                }

                * V18 FIX 2: sign-aware crossing for explain-away
                if "`criterion'" == "insignif" & `is_crossing' == 0 {
                    if `crit_s' == 0 & `prev_crit' == 0 {
                        if `sign_s' != `prev_sign' & `prev_sign' != . & `sign_s' != . {
                            local n_crossings = `n_crossings' + 2
                            local mid_s = (`prev_rd' + `rd_s') / 2
                            if `cross_closest_to_zero' == . | abs(`mid_s') < abs(`cross_closest_to_zero') {
                                local cross_closest_to_zero = `mid_s'
                                local cross_lo_at_closest    = `prev_rd'
                                local cross_hi_at_closest    = `rd_s'
                            }
                        }
                    }
                }

                if `is_crossing' == 1 {
                    local n_crossings = `n_crossings' + 1
                    if `n_crossings' == 1 {
                        local first_cross_lo = `prev_rd'
                        local first_cross_hi = `rd_s'
                    }
                    local mid_s = (`prev_rd' + `rd_s') / 2
                    if `cross_closest_to_zero' == . | abs(`mid_s') < abs(`cross_closest_to_zero') {
                        local cross_closest_to_zero = `mid_s'
                        local cross_lo_at_closest    = `prev_rd'
                        local cross_hi_at_closest    = `rd_s'
                    }
                }
            }
            local prev_rd   = `rd_s'
            local prev_crit = `crit_s'
            local prev_sign = `sign_s'
            local prev_ok   = 1
        }
        else {
            * V18 FIX 1: failed grid point breaks adjacency
            local prev_ok = 0
        }
    }

    if `n_success' == 0 {
        global sens_bstatus = "scan all failed"
        exit
    }

    *---------------------------------------------------------------------------
    * Case: no crossing anywhere
    * V18 FIX 3: check anchor before declaring robust
    *---------------------------------------------------------------------------
    if `n_crossings' == 0 {

        * If criterion is already met at rho_D = 0, threshold is zero
        if `anchor_crit' == 1 {
            global sens_rdstar  = 0
            sens_abc, rd(0) ry(`ry') rho(`rho')
            if $sens_feasible == 1 {
                sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
                sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                          atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
                capture drop _Z_hyp_
                global sens_mestar  = $sens_me
                global sens_tstar   = cond($sens_se>0, $sens_me/$sens_se, .)
            }
            global sens_bstatus = "threshold = 0"
            exit
        }

        * Genuinely robust
        local rep_rd = .
        forval s = 1/`nscan' {
            if `scan'[`s', 2] == 1 {
                local rd_s = `scan'[`s', 1]
                if `damage_dir' > 0 & (`rep_rd' == . | `rd_s' > `rep_rd') local rep_rd = `rd_s'
                if `damage_dir' < 0 & (`rep_rd' == . | `rd_s' < `rep_rd') local rep_rd = `rd_s'
            }
        }
        sens_abc, rd(`rep_rd') ry(`ry') rho(`rho')
        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        capture drop _Z_hyp_
        global sens_mestar  = $sens_me
        global sens_tstar   = cond($sens_se>0, $sens_me/$sens_se, .)
        global sens_bstatus = "robust"
        exit
    }

    *---------------------------------------------------------------------------
    * Case: exactly one crossing -> bisect within that bracket
    *---------------------------------------------------------------------------
    if `n_crossings' == 1 {
        local rd_lo = `first_cross_lo'
        local rd_hi = `first_cross_hi'

        local cross_mid = (`rd_lo' + `rd_hi') / 2
        if sign(`cross_mid') == `damage_dir' | `cross_mid' == 0 {
            global sens_bmatch = 1
        }
        else {
            global sens_bmatch = 0
        }

        sens_abc, rd(`rd_lo') ry(`ry') rho(`rho')
        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        capture drop _Z_hyp_
        local tstat_lo = $sens_me / $sens_se
        local crit_lo  = 0
        if "`criterion'" == "insignif" & abs(`tstat_lo') <  `zval'  local crit_lo = 1
        if "`criterion'" == "signflip" & sign($sens_me) != `sign0' local crit_lo = 1

        forval iter = 1/60 {
            local rd_mid = (`rd_lo' + `rd_hi') / 2
            if (`rd_hi' - `rd_lo') < `bisect_tol' continue, break

            if "`verbose'" != "" {
                di as text "  [sens_bisect bisect] model=`model' outcome=`outcome' " ///
                    "criterion=`criterion' ry=" %5.2f `ry' " iter=`iter'/60" ///
                    " rho_D_mid=" %6.3f `rd_mid'
            }

            sens_abc, rd(`rd_mid') ry(`ry') rho(`rho')
            if $sens_feasible == 0 {
                * V18 FIX 4: shrink toward the known-good (crit_lo) side
                if `crit_lo' == 0  local rd_lo = `rd_mid'
                else               local rd_hi = `rd_mid'
                continue
            }

            sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')

            timer clear 99
            timer on 99
            sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                      atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
            timer off 99
            qui timer list 99
            if "`verbose'" != "" {
                di as text "      -> fit took " %5.2f r(t99) " sec"
            }

            capture drop _Z_hyp_

            if $sens_failed == 1 | $sens_me == . {
                * V18 FIX 4: shrink toward the known-good side
                if `crit_lo' == 0  local rd_lo = `rd_mid'
                else               local rd_hi = `rd_mid'
                continue
            }

            local tstat_m = $sens_me / $sens_se
            local crit_m  = 0
            if "`criterion'" == "insignif" & abs(`tstat_m') <  `zval'   local crit_m = 1
            if "`criterion'" == "signflip" & sign($sens_me) != `sign0'  local crit_m = 1

            if `crit_m' == `crit_lo' {
                local rd_lo = `rd_mid'
            }
            else {
                local rd_hi = `rd_mid'
            }
        }

        local rd_star = (`rd_lo' + `rd_hi') / 2
        sens_abc, rd(`rd_star') ry(`ry') rho(`rho')
        if $sens_feasible == 0 {
            global sens_bstatus = "conv. infeasible"
            exit
        }
        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        capture drop _Z_hyp_

        if $sens_failed == 1 | $sens_me == . {
            global sens_bstatus = "final eval failed"
            exit
        }

        global sens_rdstar  = `rd_star'
        global sens_mestar  = $sens_me
        global sens_tstar   = $sens_me / $sens_se
        global sens_bstatus = "converged"
        exit
    }

    *---------------------------------------------------------------------------
    * Case: multiple crossings -> flag, report the threshold closest to
    * rho_D=0 (most conservative / smallest-magnitude)
    *---------------------------------------------------------------------------
    local rd_lo = `cross_lo_at_closest'
    local rd_hi = `cross_hi_at_closest'

    sens_abc, rd(`rd_lo') ry(`ry') rho(`rho')
    sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
    sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
              atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
    capture drop _Z_hyp_
    local tstat_lo = $sens_me / $sens_se
    local crit_lo  = 0
    if "`criterion'" == "insignif" & abs(`tstat_lo') <  `zval'  local crit_lo = 1
    if "`criterion'" == "signflip" & sign($sens_me) != `sign0' local crit_lo = 1

    forval iter = 1/60 {
        local rd_mid = (`rd_lo' + `rd_hi') / 2
        if (`rd_hi' - `rd_lo') < `bisect_tol' continue, break

        if "`verbose'" != "" {
            di as text "  [sens_bisect bisect-multi] model=`model' outcome=`outcome' " ///
                "criterion=`criterion' ry=" %5.2f `ry' " iter=`iter'/60" ///
                " rho_D_mid=" %6.3f `rd_mid'
        }

        sens_abc, rd(`rd_mid') ry(`ry') rho(`rho')
        if $sens_feasible == 0 {
            * V18 FIX 4: shrink toward known-good side
            if `crit_lo' == 0  local rd_lo = `rd_mid'
            else               local rd_hi = `rd_mid'
            continue
        }
        sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')

        timer clear 99
        timer on 99
        sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        timer off 99
        qui timer list 99
        if "`verbose'" != "" {
            di as text "      -> fit took " %5.2f r(t99) " sec"
        }

        capture drop _Z_hyp_
        if $sens_failed == 1 | $sens_me == . {
            * V18 FIX 4: shrink toward known-good side
            if `crit_lo' == 0  local rd_lo = `rd_mid'
            else               local rd_hi = `rd_mid'
            continue
        }
        local tstat_m = $sens_me / $sens_se
        local crit_m  = 0
        if "`criterion'" == "insignif" & abs(`tstat_m') <  `zval'   local crit_m = 1
        if "`criterion'" == "signflip" & sign($sens_me) != `sign0'  local crit_m = 1
        if `crit_m' == `crit_lo' local rd_lo = `rd_mid'
        else                     local rd_hi = `rd_mid'
    }

    local rd_star = (`rd_lo' + `rd_hi') / 2
    sens_abc, rd(`rd_star') ry(`ry') rho(`rho')
    sens_build_z, dresid(`dresid') yresid(`yresid') touse(`touse')
    sens_fit, model(`model') yvar(`yvar') dvar(`dvar') touse(`touse') ///
              atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
    capture drop _Z_hyp_

    global sens_rdstar  = `rd_star'
    global sens_mestar  = $sens_me
    global sens_tstar   = cond($sens_se>0, $sens_me/$sens_se, .)
    global sens_bstatus = "non-monotone: `n_crossings' crossings"

end


*===============================================================================
* MAIN COMMAND: sensitivity
*   UNCHANGED from Version 17.
*===============================================================================
capture program drop sensitivity
program define sensitivity, rclass
    syntax varlist(min=2) [if] [in],  ///
        Model(string)                 ///
        [ rhod(numlist)               ///
          rhoy(numlist)               ///
          seed(integer 12345)         ///
          ATMeans                     ///
          BASEoutcome(integer 1)      ///
          Level(real 95)              ///
          SUPpress                    ///
          BISect                      ///
          SIGNchange                  ///
          NSCAN(integer 51)           ///
          SCANMargin(real 0.001)      ///
          VERBose                     ///
          BENchmark(varlist) ]

    marksample touse

    local vopt = cond("`verbose'" != "", "verbose", "")

    gettoken Y varlist : varlist
    gettoken D rest    : varlist
    local X "`rest'"

    local model = lower("`model'")
    if !inlist("`model'","logit","probit","cloglog","poisson", ///
                         "nbreg","ologit","mlogit") {
        di as error "Unsupported model: `model'"
        exit 198
    }
    if "`signchange'" != "" & "`bisect'" == "" {
        di as error "signchange requires bisect"
        exit 198
    }

    local rho_d "`rhod'"
    local rho_y "`rhoy'"

    global sens_X "`X'"

    local alpha = 1 - `level'/100
    local z     = invnormal(1 - `alpha'/2)
    local atm   = 0
    if "`atmeans'" != "" local atm = 1

    if "`rho_d'" == "" local rho_d "0.1 0.2 0.3"
    if "`rho_y'" == "" local rho_y "-0.3 -0.15 0 0.15 0.3"
    local n_rho_d : word count `rho_d'
    local n_rho_y : word count `rho_y'

    set seed `seed'

    if "`suppress'" == "" {
        di as text _n "{hline 70}"
        di as text "SENSITIVITY ANALYSIS: `model' (Version 18)"
        di as text "{hline 70}"
        di as text "Treatment: `D'"
        di as text "Outcome:   `Y'"
        if "`X'" != "" di as text "Controls:  `X'"
        di as text "rho_D: `rho_d'"
        di as text "rho_Y: `rho_y'"
        di as text "Confidence level: `level'%"
    }

    *===========================================================================
    * STEP 1: outcome residuals from short model
    *===========================================================================

    if "`model'" == "logit" {
        qui logit `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', pr
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu'*(1-`mu')) if `touse'
        qui replace `Yr' = 0 if `Yr' == . & `touse'  // mu≈0 or mu≈1: treat residual as 0
    }
    else if "`model'" == "probit" {
        qui probit `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', pr
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu'*(1-`mu')) if `touse'
        qui replace `Yr' = 0 if `Yr' == . & `touse'
    }
    else if "`model'" == "cloglog" {
        qui cloglog `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', pr
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu'*(1-`mu')) if `touse'
        qui replace `Yr' = 0 if `Yr' == . & `touse'
    }
    else if "`model'" == "poisson" {
        qui poisson `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', n
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu') if `touse'
        qui replace `Yr' = 0 if `Yr' == . & `touse'  // mu≈0 edge case
    }
    else if "`model'" == "nbreg" {
        qui nbreg `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', n
        local anb = e(alpha)
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu' + `anb'*`mu'^2) if `touse'
    }
    else if "`model'" == "ologit" {
        qui ologit `Y' `D' `X' if `touse'
        qui levelsof `Y' if `touse', local(ocats)
        local usecat : word 2 of `ocats'
        if "`usecat'" == "" local usecat : word 1 of `ocats'
        tempvar pr Yr
        qui predict double `pr' if `touse', outcome(`usecat')
        qui gen double `Yr' = (`Y' == `usecat') - `pr' if `touse'
    }
    else if "`model'" == "mlogit" {
        qui mlogit `Y' `D' `X' if `touse', baseoutcome(`baseoutcome')
        qui levelsof `Y' if `touse', local(outcome_list)
        local n_outcomes : word count `outcome_list'
        foreach o of local outcome_list {
            tempvar pr_`o'
            qui predict double `pr_`o'' if `touse', outcome(`o')
            if `atm' == 1 qui margins, dydx(`D') atmeans predict(outcome(`o'))
            else          qui margins, dydx(`D') predict(outcome(`o'))
            local me_short_`o' = r(table)[1,1]
            local se_short_`o' = r(table)[2,1]
        }
    }

    if "`model'" != "mlogit" {
        qui sum `Yr' if `touse'
        if r(sd) < 1e-8 {
            di as error "Outcome residual SD near zero — check model specification"
            exit 498
        }
        qui replace `Yr' = (`Yr' - r(mean)) / r(sd) if `touse'
    }

    *===========================================================================
    * STEP 2: D residual (standardised)
    *===========================================================================
    tempvar Dr
    qui reg `D' `X' if `touse'
    qui predict double `Dr' if `touse', resid
    qui sum `Dr' if `touse'
    if r(sd) < 1e-8 {
        di as error "Treatment residual SD near zero — check that D varies and X are not collinear with D"
        exit 498
    }
    qui replace `Dr' = `Dr' / r(sd) if `touse'

    *===========================================================================
    * STEP 3: short-model AME and residual correlation
    *===========================================================================
    if "`model'" != "mlogit" {
        if      "`model'" == "logit"   qui logit   `Y' `D' `X' if `touse'
        else if "`model'" == "probit"  qui probit  `Y' `D' `X' if `touse'
        else if "`model'" == "cloglog" qui cloglog `Y' `D' `X' if `touse'
        else if "`model'" == "poisson" qui poisson `Y' `D' `X' if `touse'
        else if "`model'" == "nbreg"   qui nbreg   `Y' `D' `X' if `touse'
        else if "`model'" == "ologit"  qui ologit  `Y' `D' `X' if `touse'

        if `atm' == 1 qui margins, dydx(`D') atmeans
        else          qui margins, dydx(`D')
        local me_short = r(table)[1,1]
        local se_short = r(table)[2,1]

        qui corr `Dr' `Yr' if `touse'
        local r = r(rho)

        if "`suppress'" == "" {
            di as text _n "Short model ME: " %9.4f `me_short' ///
                " (SE: " %7.4f `se_short' ")"
            di as text "r = cor(D_resid, Y_resid) = " %7.4f `r'
        }
        * Warn if r is large — likely indicates model misspecification
        if abs(`r') > 0.1 {
            di as error _n "WARNING: r = " %6.4f `r' " (|r| > 0.1)"
            di as error "  Large r suggests the model family may not match"
            di as error "  the outcome distribution (e.g. logit on a count"
            di as error "  outcome). Check model specification before"
            di as error "  interpreting sensitivity results."
        }
    }
    else {
        foreach o of local outcome_list {
            tempvar Yr_`o'
            qui gen double `Yr_`o'' = (`Y' == `o') - `pr_`o'' if `touse'
            qui sum `Yr_`o'' if `touse'
            qui replace `Yr_`o'' = (`Yr_`o'' - r(mean)) / r(sd) if `touse'
            qui corr `Dr' `Yr_`o'' if `touse'
            local r_`o' = r(rho)
        }
        if "`suppress'" == "" {
            di as text _n "Short model MEs:"
            foreach o of local outcome_list {
                di as text "  Outcome `o': ME=" %8.4f `me_short_`o'' ///
                    "  r=" %6.4f `r_`o''
            }
        }
        * Warn if any r_o is large
        local r_warn = 0
        foreach o of local outcome_list {
            if abs(`r_`o'') > 0.1  local r_warn = 1
        }
        if `r_warn' == 1 {
            di as error _n "WARNING: one or more |r| > 0.1"
            di as error "  Large r suggests the model family may not match"
            di as error "  the outcome distribution. Check model specification"
            di as error "  before interpreting sensitivity results."
        }
    }

    *===========================================================================
    * STEP 3b: Generate fixed epsilon (used by all Z constructions)
    * Draw once and store in $sens_eps. sens_build_z reads this global
    * automatically — no need to pass it as an option through every call.
    *===========================================================================
    tempvar eps_fixed
    qui gen double `eps_fixed' = rnormal() if `touse'
    global sens_eps "`eps_fixed'"

    *===========================================================================
    * STEP 4: grid search
    *===========================================================================
    local total_cells = `n_rho_d' * `n_rho_y'

    if "`model'" != "mlogit" {

        tempname results
        matrix `results' = J(`total_cells', 9, .)
        matrix colnames `results' = rho_d rho_y ME SE CI_lo CI_hi feasible r bias

        local me_min = .
        local me_max = .

        if "`suppress'" == "" {
            di as text _n %8s "rho_D" %8s "rho_Y" %10s "ME" ///
                %8s "SE" %22s "[`level'% CI]"
            di as text "{hline 58}"
        }

        local row = 1
        foreach rd of numlist `rho_d' {
            foreach ry of numlist `rho_y' {

                sens_abc, rd(`rd') ry(`ry') rho(`r')

                if $sens_feasible == 0 {
                    mat `results'[`row',1] = `rd'
                    mat `results'[`row',2] = `ry'
                    mat `results'[`row',7] = 0
                    mat `results'[`row',8] = `r'
                    if "`suppress'" == "" {
                        di as result %8.2f `rd' %8.2f `ry' ///
                            %10s "." %8s "." %22s "infeasible"
                    }
                    local row = `row' + 1
                    continue
                }

                sens_build_z, dresid(`Dr') yresid(`Yr') touse(`touse')
                sens_fit, model(`model') yvar(`Y') dvar(`D') touse(`touse') ///
                          atm(`atm') baseoutcome(`baseoutcome')
                capture drop _Z_hyp_

                if $sens_failed == 0 & $sens_me != . {
                    local ma = $sens_me
                    local sa = $sens_se
                    local lo = `ma' - `z'*`sa'
                    local hi = `ma' + `z'*`sa'
                    mat `results'[`row',1] = `rd'
                    mat `results'[`row',2] = `ry'
                    mat `results'[`row',3] = `ma'
                    mat `results'[`row',4] = `sa'
                    mat `results'[`row',5] = `lo'
                    mat `results'[`row',6] = `hi'
                    mat `results'[`row',7]  = 1
                    mat `results'[`row',8]  = `r'
                    mat `results'[`row',9]  = `me_short' - `ma'
                    if `me_min' == . | `ma' < `me_min' local me_min = `ma'
                    if `me_max' == . | `ma' > `me_max' local me_max = `ma'
                    if "`suppress'" == "" {
                        di as result %8.2f `rd' %8.2f `ry' ///
                            %10.4f `ma' %8.4f `sa' ///
                            "  [" %7.4f `lo' "," %7.4f `hi' "]"
                    }
                }
                else {
                    mat `results'[`row',1] = `rd'
                    mat `results'[`row',2] = `ry'
                    mat `results'[`row',7] = 0
                    mat `results'[`row',8] = `r'
                    if "`suppress'" == "" {
                        di as result %8.2f `rd' %8.2f `ry' ///
                            %10s "." %8s "." %22s "failed"
                    }
                }
                local row = `row' + 1
            }
        }

        if "`suppress'" == "" {
            di as text "{hline 58}"
            di as text _n "SENSITIVITY INTERVAL"
            di as text "  [" %8.4f `me_min' ",  " %8.4f `me_max' "]"
        }

        return scalar me_short    = `me_short'
        return scalar se_short    = `se_short'
        return scalar me_min      = `me_min'
        return scalar me_max      = `me_max'
        return scalar r           = `r'
        return matrix sensitivity = `results'

    }
    else {

        foreach o of local outcome_list {
            tempname res_`o'
            matrix `res_`o'' = J(`total_cells', 8, .)
            matrix colnames `res_`o'' = rho_d rho_y ME SE CI_lo CI_hi feasible r
            local me_min_`o' = .
            local me_max_`o' = .

            if "`suppress'" == "" {
                di as text _n "{hline 60}"
                di as text "Outcome `o'"
                di as text "{hline 60}"
                di as text %8s "rho_D" %8s "rho_Y" %10s "ME" ///
                    %8s "SE" %22s "[`level'% CI]"
                di as text "{hline 58}"
            }

            local row = 1
            foreach rd of numlist `rho_d' {
                foreach ry of numlist `rho_y' {

                    sens_abc, rd(`rd') ry(`ry') rho(`r_`o'')

                    if $sens_feasible == 0 {
                        mat `res_`o''[`row',1] = `rd'
                        mat `res_`o''[`row',2] = `ry'
                        mat `res_`o''[`row',7] = 0
                        mat `res_`o''[`row',8] = `r_`o''
                        if "`suppress'" == "" {
                            di as result %8.2f `rd' %8.2f `ry' ///
                                %10s "." %8s "." %22s "infeasible"
                        }
                        local row = `row' + 1
                        continue
                    }

                    sens_build_z, dresid(`Dr') yresid(`Yr_`o'') touse(`touse')
                    sens_fit, model(`model') yvar(`Y') dvar(`D') touse(`touse') ///
                              atm(`atm') baseoutcome(`baseoutcome') outcome("`o'")
                    capture drop _Z_hyp_

                    if $sens_failed == 0 & $sens_me != . {
                        local ma = $sens_me
                        local sa = $sens_se
                        local lo = `ma' - `z'*`sa'
                        local hi = `ma' + `z'*`sa'
                        mat `res_`o''[`row',1] = `rd'
                        mat `res_`o''[`row',2] = `ry'
                        mat `res_`o''[`row',3] = `ma'
                        mat `res_`o''[`row',4] = `sa'
                        mat `res_`o''[`row',5] = `lo'
                        mat `res_`o''[`row',6] = `hi'
                        mat `res_`o''[`row',7] = 1
                        mat `res_`o''[`row',8] = `r_`o''
                        if `me_min_`o'' == . | `ma' < `me_min_`o'' local me_min_`o' = `ma'
                        if `me_max_`o'' == . | `ma' > `me_max_`o'' local me_max_`o' = `ma'
                        if "`suppress'" == "" {
                            di as result %8.2f `rd' %8.2f `ry' ///
                                %10.4f `ma' %8.4f `sa' ///
                                "  [" %7.4f `lo' "," %7.4f `hi' "]"
                        }
                    }
                    else {
                        mat `res_`o''[`row',1] = `rd'
                        mat `res_`o''[`row',2] = `ry'
                        mat `res_`o''[`row',7] = 0
                        mat `res_`o''[`row',8] = `r_`o''
                        if "`suppress'" == "" {
                            di as result %8.2f `rd' %8.2f `ry' ///
                                %10s "." %8s "." %22s "failed"
                        }
                    }
                    local row = `row' + 1
                }
            }

            if "`suppress'" == "" {
                di as text "{hline 58}"
                di as text "Sensitivity interval outcome `o': [" ///
                    %8.4f `me_min_`o'' ",  " %8.4f `me_max_`o'' "]"
            }

            return scalar me_short_`o' = `me_short_`o''
            return scalar se_short_`o' = `se_short_`o''
            return scalar me_min_`o'   = `me_min_`o''
            return scalar me_max_`o'   = `me_max_`o''
            return scalar r_`o'        = `r_`o''
            return matrix sensitivity_`o' = `res_`o''
        }

        return local outcomes = "`outcome_list'"
    }

    return local model  = "`model'"
    return local rho_d  = "`rho_d'"
    return local rho_y  = "`rho_y'"
    return scalar level = `level'

    *===========================================================================
    * STEP 5: bisection (optional)
    *===========================================================================

    if "`bisect'" == "" {
        macro drop sens_X sens_eps sens_a sens_b sens_c sens_feasible ///
                   sens_me sens_se sens_failed sens_delta      ///
                   sens_rdstar sens_mestar sens_tstar sens_bstatus
        exit
    }

    * ---------- 5a: non-mlogit ----------
    if "`model'" != "mlogit" {

        * Always run both bisections. Sign-change runs first because rho_D+
        * (the zero crossing) serves as the upper bracket for explain-away.
        * Within [0, rho_D+] the |t| function is strictly monotone: it
        * starts at t_short and decreases to 0 at the zero crossing.
        * The explain-away bisection on this sub-interval is clean and finds
        * exactly where |t| = z_alpha with ME same sign as the short model.

        tempname b_ea b_sc
        matrix `b_ea' = J(`n_rho_y', 4, .)
        matrix `b_sc' = J(`n_rho_y', 4, .)
        matrix colnames `b_ea' = rho_y rho_d_star ME_star tstat_star
        matrix colnames `b_sc' = rho_y rho_d_plus ME_plus tstat_plus

        * ---- Pass 1: Sign-change bisection ----
        * Runs first to provide the upper bracket (rho_D+) for explain-away.
        * Within [0, rho_D+] the ME stays the same sign as the short model,
        * so |t| is strictly decreasing — explain-away bisection is clean.
        * Conceptually the two criteria are independent questions; sign-change
        * runs first only for numerical bracketing, not logical dependence.
        local brow = 1
        foreach ry of numlist `rho_y' {
            sens_bisect, ry(`ry') rho(`r') meshort(`me_short')           ///
                model(`model') yvar(`Y') dvar(`D') touse(`touse')         ///
                dresid(`Dr') yresid(`Yr') zval(`z')                       ///
                criterion("signflip") atm(`atm') baseoutcome(`baseoutcome') ///
                nscan(`nscan') scanmargin(`scanmargin') `vopt'
            local zc_`brow'  = $sens_rdstar
            local zcs_`brow' = "$sens_bstatus"
            mat `b_sc'[`brow',1] = `ry'
            cap mat `b_sc'[`brow',2] = $sens_rdstar
            cap mat `b_sc'[`brow',3] = $sens_mestar
            cap mat `b_sc'[`brow',4] = $sens_tstar
            local brow = `brow' + 1
        }

        * ---- Pass 2: Explain-away bisection, capped at rho_D+ ----
        * rho_D+ serves as the upper bracket. Within [0, rho_D+] the
        * t-stat is monotone decreasing from t_short to 0, so bisection
        * finds the unique crossing where |t| = z_alpha.
        * If sign-change was robust (rho_D+ = .), the bracket is uncapped.
        di as text _n "{hline 80}"
        di as text "BISECTION (explain away): smallest rho_d* that drives AME to insignificance"
        di as text "Criterion: |ME/SE| < " %5.3f `z'
        di as text "{hline 80}"
        di as text %10s "rho_Y" %14s "rho_d*" %14s "ME at rho_d*" ///
                        %10s "t-stat" %12s "status"
        di as text "{hline 80}"

        local brow = 1
        foreach ry of numlist `rho_y' {
            local this_zc = `zc_`brow''
            local this_rb = 99
            if "`this_zc'" != "." local this_rb = `this_zc'
            sens_bisect, ry(`ry') rho(`r') meshort(`me_short')           ///
                model(`model') yvar(`Y') dvar(`D') touse(`touse')         ///
                dresid(`Dr') yresid(`Yr') zval(`z')                       ///
                criterion("insignif") atm(`atm') baseoutcome(`baseoutcome') ///
                rdbound(`this_rb') nscan(`nscan') scanmargin(`scanmargin') `vopt'

            di as result %10.4f `ry'        %14.4f $sens_rdstar  ///
                         %14.4f $sens_mestar %10.4f $sens_tstar   ///
                         %12s   "$sens_bstatus"

            mat `b_ea'[`brow',1] = `ry'
            cap mat `b_ea'[`brow',2] = $sens_rdstar
            cap mat `b_ea'[`brow',3] = $sens_mestar
            cap mat `b_ea'[`brow',4] = $sens_tstar
            local brow = `brow' + 1
        }
        di as text "{hline 80}"
        di as text "'robust' = effect remains significant across entire feasible range"
        return matrix bisect = `b_ea'

        * ---- Print sign-change results (if requested) ----
        if "`signchange'" != "" {
            di as text _n "{hline 80}"
            di as text "BISECTION (sign change, Masten & Poirier 2025):"
            di as text "  Smallest rho_d+ such that sign(AME) reverses"
            local sign_short = cond(`me_short' >= 0, "+1", "-1")
            di as text "  Short-model AME = " %9.4f `me_short' ///
                "  (sign " "`sign_short'" ")"
            di as text "{hline 80}"
            di as text %10s "rho_Y" %14s "rho_d+" %14s "ME at rho_d+" ///
                            %10s "t-stat" %12s "status"
            di as text "{hline 80}"
            local brow = 1
            foreach ry of numlist `rho_y' {
                di as result %10.4f `b_sc'[`brow',1] %14.4f `b_sc'[`brow',2] ///
                             %14.4f `b_sc'[`brow',3] %10.4f `b_sc'[`brow',4] ///
                             %12s   "`zcs_`brow''"
                local brow = `brow' + 1
            }
            di as text "{hline 80}"
            di as text "'robust' = sign never flips across entire feasible range"
            di as text _n "NOTE (Masten & Poirier 2025): sign-change point may be"
            di as text "smaller than explain-away point. Report both."
            return matrix bisect_sign = `b_sc'
        }

    }
    * ---------- 5b: Bisection for mlogit ----------
    * Same two-pass architecture as non-mlogit, repeated per outcome category.
    else {

        local n_brows = `n_outcomes' * `n_rho_y'

        * Combined matrices across all outcomes
        tempname b_ea_all b_sc_all
        matrix `b_ea_all' = J(`n_brows', 5, .)
        matrix `b_sc_all' = J(`n_brows', 5, .)
        matrix colnames `b_ea_all' = outcome rho_y rho_d_star ME_star tstat_star
        matrix colnames `b_sc_all' = outcome rho_y rho_d_plus ME_plus tstat_plus

        local grow = 1

        foreach o of local outcome_list {

            tempname b_ea_`o' b_sc_`o'
            matrix `b_ea_`o'' = J(`n_rho_y', 4, .)
            matrix `b_sc_`o'' = J(`n_rho_y', 4, .)
            matrix colnames `b_ea_`o'' = rho_y rho_d_star ME_star tstat_star
            matrix colnames `b_sc_`o'' = rho_y rho_d_plus ME_plus tstat_plus

            * ---- Pass 1: Sign-change for this outcome ----
            local brow = 1
            foreach ry of numlist `rho_y' {
                sens_bisect, ry(`ry') rho(`r_`o'') meshort(`me_short_`o'')   ///
                    model(`model') yvar(`Y') dvar(`D') touse(`touse')          ///
                    dresid(`Dr') yresid(`Yr_`o'') zval(`z')                    ///
                    criterion("signflip") atm(`atm') baseoutcome(`baseoutcome') ///
                    outcome("`o'") nscan(`nscan') scanmargin(`scanmargin') `vopt'

                local zco_`o'_`brow' = $sens_rdstar
                local zcos_`o'_`brow' = "$sens_bstatus"
                mat `b_sc_`o''[`brow',1] = `ry'
                cap mat `b_sc_`o''[`brow',2] = $sens_rdstar
                cap mat `b_sc_`o''[`brow',3] = $sens_mestar
                cap mat `b_sc_`o''[`brow',4] = $sens_tstar
                local brow = `brow' + 1
            }

            * ---- Pass 2: Explain-away for this outcome ----
            di as text _n "{hline 80}"
            di as text "BISECTION (explain away), Outcome `o'"
            di as text "{hline 80}"
            di as text %10s "rho_Y" %14s "rho_d*" %14s "ME at rho_d*" ///
                            %10s "t-stat" %12s "status"
            di as text "{hline 80}"

            local brow = 1
            local gea = (`o' - 1)*`n_rho_y' + 1
            foreach ry of numlist `rho_y' {
                local this_zco = `zco_`o'_`brow''
                local this_rb  = 99
                if "`this_zco'" != "." local this_rb = `this_zco'

                sens_bisect, ry(`ry') rho(`r_`o'') meshort(`me_short_`o'')   ///
                    model(`model') yvar(`Y') dvar(`D') touse(`touse')          ///
                    dresid(`Dr') yresid(`Yr_`o'') zval(`z')                    ///
                    criterion("insignif") atm(`atm') baseoutcome(`baseoutcome') ///
                    outcome("`o'") rdbound(`this_rb') nscan(`nscan') scanmargin(`scanmargin') `vopt'

                di as result %10.4f `ry'        %14.4f $sens_rdstar  ///
                             %14.4f $sens_mestar %10.4f $sens_tstar   ///
                             %12s   "$sens_bstatus"

                mat `b_ea_`o''[`brow',1] = `ry'
                cap mat `b_ea_`o''[`brow',2] = $sens_rdstar
                cap mat `b_ea_`o''[`brow',3] = $sens_mestar
                cap mat `b_ea_`o''[`brow',4] = $sens_tstar
                mat `b_ea_all'[`gea',1] = `o'
                mat `b_ea_all'[`gea',2] = `ry'
                cap mat `b_ea_all'[`gea',3] = $sens_rdstar
                cap mat `b_ea_all'[`gea',4] = $sens_mestar
                cap mat `b_ea_all'[`gea',5] = $sens_tstar
                local gea  = `gea'  + 1
                local brow = `brow' + 1
            }
            di as text "{hline 80}"
            return matrix bisect_`o' = `b_ea_`o''

            * ---- Print sign-change table for this outcome (if requested) ----
            if "`signchange'" != "" {
                di as text _n "{hline 80}"
                di as text "BISECTION (sign change), Outcome `o'"
                local sign_short_o = cond(`me_short_`o'' >= 0, "+1", "-1")
                di as text "  Short-model AME = " %9.4f `me_short_`o'' ///
                    "  (sign " "`sign_short_o'" ")"
                di as text "{hline 80}"
                di as text %10s "rho_Y" %14s "rho_d+" %14s "ME at rho_d+" ///
                                %10s "t-stat" %12s "status"
                di as text "{hline 80}"
                local brow = 1
                local gsc  = (`o' - 1)*`n_rho_y' + 1
                foreach ry of numlist `rho_y' {
                    di as result %10.4f `b_sc_`o''[`brow',1] ///
                                 %14.4f `b_sc_`o''[`brow',2] ///
                                 %14.4f `b_sc_`o''[`brow',3] ///
                                 %10.4f `b_sc_`o''[`brow',4] ///
                                 %12s   "`zcos_`o'_`brow''
                    mat `b_sc_all'[`gsc',1] = `o'
                    mat `b_sc_all'[`gsc',2] = `ry'
                    cap mat `b_sc_all'[`gsc',3] = `b_sc_`o''[`brow',2]
                    cap mat `b_sc_all'[`gsc',4] = `b_sc_`o''[`brow',3]
                    cap mat `b_sc_all'[`gsc',5] = `b_sc_`o''[`brow',4]
                    local gsc  = `gsc'  + 1
                    local brow = `brow' + 1
                }
                di as text "{hline 80}"
                return matrix bisect_sign_`o' = `b_sc_`o''
            }

            local grow = `grow' + `n_rho_y'

        } // end foreach o

        return matrix bisect = `b_ea_all'
        if "`signchange'" != "" return matrix bisect_sign = `b_sc_all'

    } // end mlogit bisection

    * (globals cleaned up at end of program)

    *===========================================================================
    * STEP 6: benchmark analysis (optional)
    * For each variable in benchmark(varlist):
    *   1. Drop the variable from the short model and re-estimate.
    *   2. Compute rho_D_k = cor(Dr_k, Dr) where Dr_k is the variable
    *      residualised on the remaining controls.
    *   3. Compute rho_Y_k = cor(Dr_k, Yr) — partial correlation of the
    *      benchmarking variable with the outcome residual.
    *   4. Evaluate the AME at (rho_D_k, rho_Y_k) by constructing Z
    *      and fitting the long model.
    *   5. Report rho_D_k, rho_Y_k, the adjusted AME, and the implied
    *      bias — how much the short-model AME would shift if a confounder
    *      as strong as this variable were omitted.
    *===========================================================================

    if "`benchmark'" == "" {
        macro drop sens_X sens_eps sens_a sens_b sens_c sens_feasible ///
                   sens_me sens_se sens_failed sens_delta             ///
                   sens_rdstar sens_mestar sens_tstar sens_bstatus
        exit
    }

    if "`model'" != "mlogit" {

        di as text _n "{hline 80}"
        di as text "BENCHMARK ANALYSIS: sensitivity at observed-covariate strength"
        di as text "{hline 80}"
        di as text %14s "variable" %10s "rho_D" %10s "rho_Y" ///
                        %12s "ME" %10s "SE" %10s "t-stat" ///
                        %10s "bias" %10s "feasible"
        di as text "{hline 80}"

        local n_bench : word count `benchmark'
        tempname bench_mat
        matrix `bench_mat' = J(`n_bench', 7, .)
        matrix colnames `bench_mat' = rho_D rho_Y ME SE tstat bias feasible

        local bk = 1
        foreach bvar of varlist `benchmark' {

            local X_minus ""
            foreach xv of local X {
                if "`xv'" != "`bvar'" {
                    local X_minus "`X_minus' `xv'"
                }
            }

            * Residualise benchmark variable on remaining controls
            tempvar W_r Dr_b Yr_b
            if "`X_minus'" != "" {
                qui reg `bvar' `X_minus' if `touse'
                qui predict double `W_r' if `touse', resid
            }
            else {
                qui gen double `W_r' = `bvar' if `touse'
            }
            qui sum `W_r' if `touse'
            if r(sd) < 1e-10 {
                di as result %14s "`bvar'" "  (zero variance — skipped)"
                local bk = `bk' + 1
                continue
            }
            qui replace `W_r' = `W_r' / r(sd) if `touse'

            * Treatment residual from model excluding benchmark variable
            qui reg `D' `X_minus' if `touse'
            qui predict double `Dr_b' if `touse', resid
            qui sum `Dr_b' if `touse'
            if r(sd) > 1e-10 {
                qui replace `Dr_b' = `Dr_b' / r(sd) if `touse'
            }

            /** Outcome residual from OLS excluding benchmark variable
            * Using OLS gives non-trivial rho_Y (unlike logit Pearson residuals
            * which are near-zero by score orthogonality in the main model)
            qui reg `Y' `D' `X_minus' if `touse'
            qui predict double `Yr_b' if `touse', resid
            qui sum `Yr_b' if `touse'
            if r(sd) > 1e-10 {
                qui replace `Yr_b' = `Yr_b' / r(sd) if `touse'
            }*/
			
			* Outcome residual from reduced nonlinear model excluding benchmark variable
            tempvar mu_b Yr_b
            if "`model'" == "logit" {
                qui logit `Y' `D' `X_minus' if `touse'
                qui predict double `mu_b' if `touse', pr
                qui gen double `Yr_b' = (`Y' - `mu_b') / sqrt(`mu_b'*(1-`mu_b')) if `touse'
                qui replace `Yr_b' = 0 if `Yr_b' == . & `touse'
            }
            else if "`model'" == "probit" {
                qui probit `Y' `D' `X_minus' if `touse'
                qui predict double `mu_b' if `touse', pr
                qui gen double `Yr_b' = (`Y' - `mu_b') / sqrt(`mu_b'*(1-`mu_b')) if `touse'
                qui replace `Yr_b' = 0 if `Yr_b' == . & `touse'
            }
            else if "`model'" == "cloglog" {
                qui cloglog `Y' `D' `X_minus' if `touse'
                qui predict double `mu_b' if `touse', pr
                qui gen double `Yr_b' = (`Y' - `mu_b') / sqrt(`mu_b'*(1-`mu_b')) if `touse'
                qui replace `Yr_b' = 0 if `Yr_b' == . & `touse'
            }
            else if "`model'" == "poisson" {
                qui poisson `Y' `D' `X_minus' if `touse'
                qui predict double `mu_b' if `touse', n
                qui gen double `Yr_b' = (`Y' - `mu_b') / sqrt(`mu_b') if `touse'
                qui replace `Yr_b' = 0 if `Yr_b' == . & `touse'
            }
            else if "`model'" == "nbreg" {
                qui nbreg `Y' `D' `X_minus' if `touse'
                qui predict double `mu_b' if `touse', n
                local anb_b = e(alpha)
                qui gen double `Yr_b' = (`Y' - `mu_b') / ///
                    sqrt(`mu_b' + `anb_b'*`mu_b'^2) if `touse'
            }
            else if "`model'" == "ologit" {
                qui ologit `Y' `D' `X_minus' if `touse'
                qui levelsof `Y' if `touse', local(ocats_b)
                local usecat_b : word 2 of `ocats_b'
                if "`usecat_b'" == "" local usecat_b : word 1 of `ocats_b'
                qui predict double `mu_b' if `touse', outcome(`usecat_b')
                qui gen double `Yr_b' = (`Y' == `usecat_b') - `mu_b' if `touse'
            }
            qui sum `Yr_b' if `touse'
            if r(sd) < 1e-10 {
                di as result %14s "`bvar'" "  (zero residual SD — skipped)"
                local bk = `bk' + 1
                continue
            }
            qui replace `Yr_b' = (`Yr_b' - r(mean)) / r(sd) if `touse'

            * Benchmark partial correlations
            qui corr `W_r' `Dr_b' if `touse'
            local rho_D_k = r(rho)
            qui corr `W_r' `Yr_b' if `touse'
            local rho_Y_k = r(rho)

            local rho_D_k = max(min(`rho_D_k', 0.999), -0.999)
            local rho_Y_k = max(min(`rho_Y_k', 0.999), -0.999)

            * Check feasibility and evaluate adjusted AME
            sens_abc, rd(`rho_D_k') ry(`rho_Y_k') rho(`r')

            if $sens_feasible == 0 {
                di as result %14s "`bvar'" ///
                    %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                    %12s "." %10s "." %10s "." %10s "." %10s "infeasible"
                matrix `bench_mat'[`bk',1] = `rho_D_k'
                matrix `bench_mat'[`bk',2] = `rho_Y_k'
                matrix `bench_mat'[`bk',7] = 0
                local bk = `bk' + 1
                continue
            }

            sens_build_z, dresid(`Dr') yresid(`Yr') touse(`touse')
            sens_fit, model(`model') yvar(`Y') dvar(`D') ///
                      touse(`touse') atm(`atm') baseoutcome(`baseoutcome')
            capture drop _Z_hyp_

            if $sens_failed == 0 & $sens_me != . {
                local bk_me   = $sens_me
                local bk_se   = $sens_se
                local bk_t    = `bk_me' / `bk_se'
                local bk_bias = `me_short' - `bk_me'
                di as result %14s "`bvar'" ///
                    %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                    %12.4f `bk_me' %10.4f `bk_se' %10.4f `bk_t' ///
                    %10.4f `bk_bias' %10s "feasible"
                matrix `bench_mat'[`bk',1] = `rho_D_k'
                matrix `bench_mat'[`bk',2] = `rho_Y_k'
                matrix `bench_mat'[`bk',3] = `bk_me'
                matrix `bench_mat'[`bk',4] = `bk_se'
                matrix `bench_mat'[`bk',5] = `bk_t'
                matrix `bench_mat'[`bk',6] = `bk_bias'
                matrix `bench_mat'[`bk',7] = 1
            }
            else {
                di as result %14s "`bvar'" ///
                    %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                    %12s "." %10s "." %10s "." %10s "." %10s "failed"
                matrix `bench_mat'[`bk',1] = `rho_D_k'
                matrix `bench_mat'[`bk',2] = `rho_Y_k'
                matrix `bench_mat'[`bk',7] = 0
            }
            local bk = `bk' + 1
        }
        di as text "{hline 80}"
        di as text "Interpretation: 'bias' = ME_short - ME_adjusted;"
        di as text "  if a confounder as strong as this variable were omitted,"
        di as text "  the short-model AME would shift by this amount."
        matrix rownames `bench_mat' = `benchmark'
        return matrix benchmark = `bench_mat'
    }
    else {

        foreach o of local outcome_list {

            di as text _n "{hline 80}"
            di as text "BENCHMARK ANALYSIS: Outcome `o'"
            di as text "{hline 80}"
            di as text %14s "variable" %10s "rho_D" %10s "rho_Y" ///
                            %12s "ME" %10s "SE" %10s "t-stat" ///
                            %10s "bias" %10s "feasible"
            di as text "{hline 80}"

            local n_bench : word count `benchmark'
            tempname bench_mat_`o'
            matrix `bench_mat_`o'' = J(`n_bench', 7, .)
            matrix colnames `bench_mat_`o'' = rho_D rho_Y ME SE tstat bias feasible

            local bk = 1
            foreach bvar of varlist `benchmark' {

                local X_minus ""
                foreach xv of local X {
                    if "`xv'" != "`bvar'" local X_minus "`X_minus' `xv'"
                }

                tempvar W_r Dr_b pr_b_o Yr_b
                if "`X_minus'" != "" {
                    qui reg `bvar' `X_minus' if `touse'
                    qui predict double `W_r' if `touse', resid
                }
                else {
                    qui gen double `W_r' = `bvar' if `touse'
                }
                qui sum `W_r' if `touse'
                if r(sd) < 1e-10 {
                    di as result %14s "`bvar'" "  (zero variance — skipped)"
                    local bk = `bk' + 1
                    continue
                }
                qui replace `W_r' = `W_r' / r(sd) if `touse'

                qui reg `D' `X_minus' if `touse'
                qui predict double `Dr_b' if `touse', resid
                qui sum `Dr_b' if `touse'
                if r(sd) > 1e-10 {
                    qui replace `Dr_b' = `Dr_b' / r(sd) if `touse'
                }

                * Re-estimate mlogit without benchmark variable
                qui mlogit `Y' `D' `X_minus' if `touse', baseoutcome(`baseoutcome')
                qui predict double `pr_b_o' if `touse', outcome(`o')
                qui gen double `Yr_b' = (`Y' == `o') - `pr_b_o' if `touse'
                qui sum `Yr_b' if `touse'
                if r(sd) < 1e-10 {
                    di as result %14s "`bvar'" "  (zero residual SD — skipped)"
                    local bk = `bk' + 1
                    continue
                }
                qui replace `Yr_b' = (`Yr_b' - r(mean)) / r(sd) if `touse'

                qui corr `W_r' `Dr_b' if `touse'
                local rho_D_k = r(rho)
                qui corr `W_r' `Yr_b' if `touse'
                local rho_Y_k = r(rho)
                local rho_D_k = max(min(`rho_D_k',  0.999), -0.999)
                local rho_Y_k = max(min(`rho_Y_k',  0.999), -0.999)

                * V18 FIX 5: use MAIN-MODEL category-specific r_o, not
                * reduced-model r_bench. The benchmark correlations come
                * from the reduced model, but Z is built from the main-
                * model residuals, so the loading formula must use the r
                * that matches those residuals.
                sens_abc, rd(`rho_D_k') ry(`rho_Y_k') rho(`r_`o'')

                if $sens_feasible == 0 {
                    di as result %14s "`bvar'" ///
                        %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                        %12s "." %10s "." %10s "." %10s "." %10s "infeasible"
                    matrix `bench_mat_`o''[`bk',1] = `rho_D_k'
                    matrix `bench_mat_`o''[`bk',2] = `rho_Y_k'
                    matrix `bench_mat_`o''[`bk',7] = 0
                    local bk = `bk' + 1
                    continue
                }

                sens_build_z, dresid(`Dr') yresid(`Yr_`o'') touse(`touse')
                sens_fit, model(`model') yvar(`Y') dvar(`D') ///
                          touse(`touse') atm(`atm')           ///
                          baseoutcome(`baseoutcome') outcome("`o'")
                capture drop _Z_hyp_

                if $sens_failed == 0 & $sens_me != . {
                    local bk_me   = $sens_me
                    local bk_se   = $sens_se
                    local bk_t    = `bk_me' / `bk_se'
                    local bk_bias = `me_short_`o'' - `bk_me'
                    di as result %14s "`bvar'" ///
                        %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                        %12.4f `bk_me' %10.4f `bk_se' %10.4f `bk_t' ///
                        %10.4f `bk_bias' %10s "feasible"
                    matrix `bench_mat_`o''[`bk',1] = `rho_D_k'
                    matrix `bench_mat_`o''[`bk',2] = `rho_Y_k'
                    matrix `bench_mat_`o''[`bk',3] = `bk_me'
                    matrix `bench_mat_`o''[`bk',4] = `bk_se'
                    matrix `bench_mat_`o''[`bk',5] = `bk_t'
                    matrix `bench_mat_`o''[`bk',6] = `bk_bias'
                    matrix `bench_mat_`o''[`bk',7] = 1
                }
                else {
                    di as result %14s "`bvar'" ///
                        %10.4f `rho_D_k' %10.4f `rho_Y_k' ///
                        %12s "." %10s "." %10s "." %10s "." %10s "failed"
                    matrix `bench_mat_`o''[`bk',1] = `rho_D_k'
                    matrix `bench_mat_`o''[`bk',2] = `rho_Y_k'
                    matrix `bench_mat_`o''[`bk',7] = 0
                }
                local bk = `bk' + 1
            }
            di as text "{hline 80}"
            di as text "Interpretation: 'bias' = ME_short - ME_adjusted;"
            di as text "  if a confounder as strong as this variable were omitted,"
            di as text "  the short-model AME would shift by this amount."
            matrix rownames `bench_mat_`o'' = `benchmark'
            return matrix benchmark_`o' = `bench_mat_`o''
        }
    }

    macro drop sens_X sens_eps sens_a sens_b sens_c sens_feasible ///
               sens_me sens_se sens_failed sens_delta             ///
               sens_rdstar sens_mestar sens_tstar sens_bstatus

end


log using "/Users/user/Library/CloudStorage/Dropbox/Econometrics/politicalBehaviour/simulations.smcl", replace


*===============================================================================
* sensitivity_example
*   UNCHANGED from Version 17 (except header text noting the new bisect
*   engine), original 3-model illustration (logit, Poisson, mlogit).
*===============================================================================
capture program drop sensitivity_example
program define sensitivity_example

    di as text _n "{hline 70}"
    di as text "SENSITIVITY ANALYSIS EXAMPLES (Version 18)"
    di as text "{hline 70}"

    clear
    set seed 54321
    set obs 20000

    gen X1 = rnormal()
    gen X2 = rnormal()
    gen D  = (0.3*X1 + rnormal() > 0)

    qui reg D X1 X2
    qui predict Dr, resid
    qui sum Dr
    qui replace Dr = Dr/r(sd)
    gen Z = 0.3*Dr + sqrt(1-0.09)*rnormal()
    drop Dr

    gen eta_b = -0.5 + 0.5*D + 0.3*X1 + 0.2*X2 + 0.5*Z
    gen Y_bin = runiform() < invlogit(eta_b)

    gen eta_c = 0.5 + 0.3*D + 0.2*X1 + 0.1*X2 + 0.3*Z
    gen Y_cnt = rpoisson(exp(eta_c))

    gen eta2  = -0.5 + 0.5*D + 0.2*X1 + 0.3*X2 + 0.4*Z
    gen eta3  = -1.0 + 0.7*D + 0.1*X1 + 0.4*X2 + 0.5*Z
    gen denom = 1 + exp(eta2) + exp(eta3)
    gen u     = runiform()
    gen Y_mul = cond(u < 1/denom, 1, cond(u < (1+exp(eta2))/denom, 2, 3))

    di as text _n "{hline 70}"
    di as text "Example 1: Binary Logit"
    di as text "{hline 70}"
    sensitivity Y_bin D X1 X2, model(logit) ///
        rhod(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        rhoy(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        bisect signchange benchmark(X1 X2)

    di as text _n "{hline 70}"
    di as text "Example 2: Poisson"
    di as text "{hline 70}"
    sensitivity Y_cnt D X1 X2, model(poisson) ///
        rhod(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        rhoy(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        bisect signchange benchmark(X1 X2)

    di as text _n "{hline 70}"
    di as text "Example 3: Multinomial Logit"
    di as text "{hline 70}"
    sensitivity Y_mul D X1 X2, model(mlogit) baseoutcome(1) ///
        rhod(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        rhoy(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5) ///
        bisect signchange benchmark(X1 X2)

end

sensitivity_example



*===============================================================================
* sensitivity_dgp_combined
*
* PURPOSE: Three-DGP validation. For each DGP, data are generated ONCE and
*   the sensitivity command is called TWICE on the same dataset:
*     Pass 1: benchmark(X1 X2)       — observed covariates only
*     Pass 2: benchmark(X1 X2 Z)     — adds true confounder for validation
*
*   This ensures the threshold table and benchmark validation table are
*   directly comparable (same data, same epsilon draw via seed()).
*
*   DGP 1 — No confounder: Z has zero effect on Y.
*            Expected: robust thresholds; Z benchmark shows near-zero rho_Y.
*   DGP 2 — Weak confounder: Z coeff=0.2, rho_D≈0.3.
*            Expected: moderate fragility; Z benchmark shows small bias.
*   DGP 3 — Strong confounder: Z coeff=0.8, rho_D≈0.6, true D coeff=0.3.
*            Expected: high fragility; Z benchmark shows large bias.
*   DGP 4 — Zero true effect: spurious significance driven entirely by
*            confounding; threshold should fall below the true rho_D.
*
*   UNCHANGED from the original Version 17 script. Calls only the public
*   `sensitivity` command, whose signature is identical in the
*   validate-then-bisect version (only sens_bisect's internals differ,
*   plus the new optional nscan()/scanmargin()/verbose options, which
*   default to behavior matching the original when omitted, as used here).
*===============================================================================
capture program drop sensitivity_dgp_combined
program define sensitivity_dgp_combined

    local rhod_grid "-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7"
    local rhoy_grid "-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5"

    *---------------------------------------------------------------------------
    * DGP 1: No confounder
    *---------------------------------------------------------------------------
    di as text _n "{hline 70}"
    di as text "DGP 1: No omitted confounder"
    di as text "  True model: Y ~ D + X1 + X2  (Z has zero effect on Y)"
    di as text "  True treatment coeff (logit): 0.5"
    di as text "  Z has rho_D approx 0.5 but zero effect on Y"
    di as text "  Expected: robust thresholds; Z benchmark near-zero rho_Y"
    di as text "{hline 70}"

    clear
    set obs 20000
    gen X1 = rnormal()
    gen X2 = rnormal()
    gen D  = (0.3*X1 + rnormal() > 0)

    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr/r(sd)
    gen Z = 0.5*Dr + sqrt(1-0.25)*rnormal()
    drop Dr

    gen eta = -0.5 + 0.5*D + 0.3*X1 + 0.2*X2
    gen Y   = runiform() < invlogit(eta)

    di as text _n "--- Pass 1: grid + bisection + benchmark(X1 X2) ---"
    sensitivity Y D X1 X2, model(logit) seed(11111) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2)

    di as text _n "--- Pass 2: benchmark(X1 X2 Z) — validation ---"
    sensitivity Y D X1 X2, model(logit) seed(11111) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2 Z) suppress

	 
    global me_short_dgp1 = r(me_short)


    *---------------------------------------------------------------------------
    * DGP 2: Weak confounder
    *---------------------------------------------------------------------------
    di as text _n "{hline 70}"
    di as text "DGP 2: Weak omitted confounder"
    di as text "  True model: Y ~ D + X1 + X2 + 0.2*Z"
    di as text "  True treatment coeff (logit): 0.5"
    di as text "  Z: rho_D approx 0.3, effect on Y = 0.2"
    di as text "  Expected: moderate fragility; Z benchmark shows small bias"
    di as text "{hline 70}"

    clear
    set obs 20000
    gen X1 = rnormal()
    gen X2 = rnormal()
    gen D  = (0.3*X1 + rnormal() > 0)

    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr/r(sd)
    gen Z = 0.3*Dr + sqrt(1-0.09)*rnormal()
    drop Dr

    gen eta = -0.5 + 0.5*D + 0.3*X1 + 0.2*X2 + 0.2*Z
    gen Y   = runiform() < invlogit(eta)

    di as text _n "--- Pass 1: grid + bisection + benchmark(X1 X2) ---"
    sensitivity Y D X1 X2, model(logit) seed(22222) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2)

    di as text _n "--- Pass 2: benchmark(X1 X2 Z) — validation ---"
    sensitivity Y D X1 X2, model(logit) seed(22222) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2 Z) suppress

	global me_short_dgp2 = r(me_short)
	
    *---------------------------------------------------------------------------
    * DGP 3: Strong confounder
    *---------------------------------------------------------------------------
    di as text _n "{hline 70}"
    di as text "DGP 3: Strong omitted confounder"
    di as text "  True model: Y ~ D + X1 + X2 + 0.8*Z"
    di as text "  True treatment coeff (logit): 0.3  (deliberately weak)"
    di as text "  Z: rho_D approx 0.6, effect on Y = 0.8"
    di as text "  Expected: high fragility; Z benchmark shows large bias"
    di as text "{hline 70}"

    clear
    set obs 20000
    gen X1 = rnormal()
    gen X2 = rnormal()
    gen D  = (0.3*X1 + rnormal() > 0)

    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr/r(sd)
    gen Z = 0.6*Dr + sqrt(1-0.36)*rnormal()
    drop Dr

    gen eta = -0.5 + 0.3*D + 0.3*X1 + 0.2*X2 + 0.8*Z
    gen Y   = runiform() < invlogit(eta)

    di as text _n "--- Pass 1: grid + bisection + benchmark(X1 X2) ---"
    sensitivity Y D X1 X2, model(logit) seed(33333) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2)

    di as text _n "--- Pass 2: benchmark(X1 X2 Z) — validation ---"
    sensitivity Y D X1 X2, model(logit) seed(33333) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2 Z) suppress

	global me_short_dgp3 = r(me_short)

	*---------------------------------------------------------------------------
    * DGP 4: Zero true effect, spurious significance from confounding
    *---------------------------------------------------------------------------
    di as text _n "{hline 70}"
    di as text "DGP 4: Zero true effect — spurious significance"
    di as text "  True model: Y ~ logit(-0.5 + 0.0*D + 0.3*X1 + 0.2*X2 + 0.8*Z)"
    di as text "  True treatment coeff (logit): 0.0  (no causal effect)"
    di as text "  Z: rho_D approx 0.4, effect on Y = 0.8"
    di as text "  Expected: significant short-model AME that vanishes under"
    di as text "  confounding adjustment; threshold below true rho_D"
    di as text "{hline 70}"

    clear
    set obs 20000
    gen X1 = rnormal()
    gen X2 = rnormal()
    gen D  = (0.3*X1 + rnormal() > 0)

    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr/r(sd)
    gen Z = 0.4*Dr + sqrt(1-0.16)*rnormal()
    drop Dr

    gen eta = -0.5 + 0.0*D + 0.3*X1 + 0.2*X2 + 0.8*Z
    gen Y   = runiform() < invlogit(eta)

    di as text _n "--- Pass 1: grid + bisection + benchmark(X1 X2) ---"
    sensitivity Y D X1 X2, model(logit) seed(44444) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2)

    di as text _n "--- Pass 2: benchmark(X1 X2 Z) — validation ---"
    sensitivity Y D X1 X2, model(logit) seed(44444) ///
        rhod(`rhod_grid') rhoy(`rhoy_grid')          ///
        bisect signchange benchmark(X1 X2 Z) suppress

	 global me_short_dgp4 = r(me_short)

end

sensitivity_dgp_combined


*===============================================================================
* true_ame
* Computes the true population AME for each DGP by Monte Carlo integration.
* Uses a very large sample (N=500,000) to approximate the expectation
* E[Lambda'(alpha + beta*D + gamma1*X1 + gamma2*X2 + delta*Z)]
* The true AME = beta * E[Lambda'(eta)] where eta is the true linear index.
*
* UNCHANGED from the original Version 17 script. Does not call the
* sensitivity package at all — pure Monte Carlo integration plus a
* readback of the me_short_dgp1..4 globals set by sensitivity_dgp_combined.
*===============================================================================

capture program drop true_ame
program define true_ame

    di as text _n "{hline 70}"
    di as text "TRUE POPULATION AME BY MONTE CARLO INTEGRATION"
    di as text "N = 500,000 draws from the true DGP"
    di as text "{hline 70}"

    *---------------------------------------------------------------------------
    * DGP 1
    *---------------------------------------------------------------------------
    di as text _n "DGP 1: No confounder"
    clear
    set obs 500000
    gen X1  = rnormal()
    gen X2  = rnormal()
    gen D   = (0.3*X1 + rnormal() > 0)
    gen eta = -0.5 + 0.5*D + 0.3*X1 + 0.2*X2
    gen lp  = invlogit(eta)
    gen lpd = lp * (1 - lp)
    qui sum lpd
    local true_ame1 = 0.5 * r(mean)
    di as result "  True AME (D) = " %8.4f `true_ame1'

    *---------------------------------------------------------------------------
    * DGP 2
    *---------------------------------------------------------------------------
    di as text _n "DGP 2: Weak confounder"
    clear
    set obs 500000
    gen X1  = rnormal()
    gen X2  = rnormal()
    gen D   = (0.3*X1 + rnormal() > 0)
    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr / r(sd)
    gen Z   = 0.3*Dr + sqrt(1-0.09)*rnormal()
    drop Dr
    gen eta = -0.5 + 0.5*D + 0.3*X1 + 0.2*X2 + 0.2*Z
    gen lp  = invlogit(eta)
    gen lpd = lp * (1 - lp)
    qui sum lpd
    local true_ame2 = 0.5 * r(mean)
    di as result "  True AME (D) = " %8.4f `true_ame2'

    *---------------------------------------------------------------------------
    * DGP 3
    *---------------------------------------------------------------------------
    di as text _n "DGP 3: Strong confounder"
    clear
    set obs 500000
    gen X1  = rnormal()
    gen X2  = rnormal()
    gen D   = (0.3*X1 + rnormal() > 0)
    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr / r(sd)
    gen Z   = 0.6*Dr + sqrt(1-0.36)*rnormal()
    drop Dr
    gen eta = -0.5 + 0.3*D + 0.3*X1 + 0.2*X2 + 0.8*Z
    gen lp  = invlogit(eta)
    gen lpd = lp * (1 - lp)
    qui sum lpd
    local true_ame3 = 0.3 * r(mean)
    di as result "  True AME (D) = " %8.4f `true_ame3'

    *---------------------------------------------------------------------------
    * DGP 4
    *---------------------------------------------------------------------------
    di as text _n "DGP 4: Zero true effect"
    clear
    set obs 500000
    gen X1  = rnormal()
    gen X2  = rnormal()
    gen D   = (0.3*X1 + rnormal() > 0)
    qui reg D X1 X2
    qui predict double Dr, resid
    qui sum Dr
    qui replace Dr = Dr / r(sd)
    gen Z   = 0.4*Dr + sqrt(1-0.16)*rnormal()
    drop Dr
    gen eta = -0.5 + 0.0*D + 0.3*X1 + 0.2*X2 + 0.8*Z
    gen lp  = invlogit(eta)
    gen lpd = lp * (1 - lp)
    qui sum lpd
    local true_ame4 = 0.0 * r(mean)   // exactly zero by construction
    di as result "  True AME (D) = " %8.4f `true_ame4'

    *---------------------------------------------------------------------------
    * Summary — short-model AMEs read from globals set by sensitivity_dgp_combined
    *---------------------------------------------------------------------------
    * Check globals are available
    foreach g in me_short_dgp1 me_short_dgp2 me_short_dgp3 me_short_dgp4 {
        if "${`g'}" == "" {
            di as error "Global `g' not found — run sensitivity_dgp_combined first"
            exit 198
        }
    }

    local sme1 = ${me_short_dgp1}
    local sme2 = ${me_short_dgp2}
    local sme3 = ${me_short_dgp3}
    local sme4 = ${me_short_dgp4}

    di as text _n "{hline 70}"
    di as text "SUMMARY"
    di as text "{hline 70}"
    di as text %25s "DGP" %12s "True AME" %18s "Short-model AME" %12s "Bias"
    di as text "{hline 67}"
    di as result %25s "1 (No confounder)"        %12.4f `true_ame1' ///
                 %18.4f `sme1' %12.4f (`sme1' - `true_ame1')
    di as result %25s "2 (Weak confounder)"      %12.4f `true_ame2' ///
                 %18.4f `sme2' %12.4f (`sme2' - `true_ame2')
    di as result %25s "3 (Strong confounder)"    %12.4f `true_ame3' ///
                 %18.4f `sme3' %12.4f (`sme3' - `true_ame3')
    di as result %25s "4 (Zero true effect)"     %12.4f `true_ame4' ///
                 %18.4f `sme4' %12.4f (`sme4' - `true_ame4')
    di as text "{hline 67}"
    di as text "Note: Short-model AMEs are globals set by sensitivity_dgp_combined."
    di as text "Bias = Short-model AME - True AME."

end

true_ame


*===============================================================================
* simulation_comparison_v3
*
* PURPOSE: Compare three sensitivity analysis methods across a grid of
*   confounding strengths. Tests the theoretical prediction that the
*   CH linear approximation error grows with delta_Z (O(delta^2) remainder)
*   while the synthetic method remains accurate throughout.
*
* NOTE ON SCOPE: this simulation does NOT call the sensitivity package
* (sensitivity / sens_bisect / sens_fit / sens_abc / sens_build_z) at
* all. Its "Synthetic" method reimplements the a/b/c synthetic-confounder
* construction inline, independently of sens_abc. It is therefore
* UNAFFECTED by the validate-then-bisect redesign of sens_bisect, and is
* included here unchanged purely so this single file reproduces every
* simulation from the original script. If you only care about results
* that depend on sens_bisect's behavior, this section can be skipped.
*
* DESIGN:
*   6 cells: 2 treatment effect levels x 3 confounding strengths.
*   Everything else fixed: intc=-0.5, scale=1, nobs=1000.
*
*   Cell 1: beta_D=0,   delta_Z=0.2, rho_D=0.3  [null, very weak]
*   Cell 2: beta_D=0,   delta_Z=0.5, rho_D=0.4  [null, moderate]
*   Cell 3: beta_D=0,   delta_Z=0.8, rho_D=0.6  [null, strong]
*   Cell 4: beta_D=0.5, delta_Z=0.2, rho_D=0.3  [positive, very weak]
*   Cell 5: beta_D=0.5, delta_Z=0.5, rho_D=0.4  [positive, moderate]
*   Cell 6: beta_D=0.5, delta_Z=0.8, rho_D=0.6  [positive, strong]
*
* KEY PREDICTION:
*   MAE advantage of synthetic over CH should grow monotonically with
*   delta_Z, consistent with the O(delta^2) remainder in the OVB lemma.
*   CH bias should be approximately zero at delta_Z=0.2 and grow
*   quadratically as delta_Z increases to 0.5 and 0.8.
*
* PRIMARY METRIC:  MAE = mean|adjusted AME - true AME|
* SECONDARY:       Bias = mean(adjusted AME - true AME)
*                  Direction error rate
*
* METHODS:
*   CH-OLS:    short OLS + OLS residuals + CH formula
*   CH-adapted: short logit AME + logit Pearson residuals + CH formula
*   Synthetic: short logit + logit Pearson residuals + exact re-estimation
*
* GROUND TRUTH: true population AME from logit on N=200,000 draws.
* NO conditioning on short-model significance.
*
* SEEDING CAVEAT (unchanged from the original — flagged, not fixed here):
*   A single `set seed 12345` is called once, immediately before
*   `sim_v3, nsim(10000) nobs(5000)` at the bottom of this section. Every
*   random draw inside sim_v3 (the population true-AME draw for each
*   cell, then that cell's nsim replications) consumes one continuous
*   stream in a fixed order. The full run is exactly reproducible, but
*   re-running a single cell in isolation, or changing nsim/nobs, will
*   NOT reproduce any individual cell's original numbers, because the
*   RNG state at the start of a later cell depends on exactly how many
*   draws every earlier cell consumed. If per-cell reproducibility is
*   wanted, each cell would need its own seed (e.g. derived from a base
*   seed plus the cell index) — not changed here to keep this section a
*   faithful, unmodified port of the original.
*===============================================================================

capture program drop ch_ols_adjusted
program define ch_ols_adjusted
    syntax , b(real) se(real) n(real) r2dz(real) r2yz(real) ///
             [ ncontrols(integer 2) ]

    global ch_t_adj = .
    global ch_b_adj = .

    if `r2dz' >= 1 | `r2yz' >= 1 | `r2dz' < 0 | `r2yz' < 0 exit
    if `se' <= 0 exit

    local dof    = `n' - `ncontrols' - 2
    if `dof' <= 0 exit

    local t_short   = `b' / `se'
    local sign_t    = sign(`t_short')
    local denom_r2  = max(1e-10, 1 - `r2dz')
    local bf        = sqrt(`r2yz' * `r2dz' / `denom_r2')
    local bias_on_t = `sign_t' * `bf' * sqrt(`dof')

    global ch_t_adj = `t_short' - `bias_on_t'
    global ch_b_adj = `b' - `sign_t' * `bf' * `se' * sqrt(`dof')

end


capture program drop compute_true_ame
program define compute_true_ame
    syntax , intc(real) beta(real) g1(real) g2(real) ///
             delta(real) rho0(real)

    preserve
    quietly {
        clear
        set obs 500000
        gen X1 = rnormal()
        gen X2 = rnormal()
        gen D  = (0.3*X1 + rnormal() > 0)

        reg D X1 X2
        predict double Dr0, resid
        sum Dr0
        replace Dr0 = Dr0 / r(sd)
        gen double Z = `rho0'*Dr0 + sqrt(1-`rho0'^2)*rnormal()
        drop Dr0

        gen double eta = `intc' + `beta'*D + `g1'*X1 + `g2'*X2 + `delta'*Z
        gen byte   Y   = runiform() < invlogit(eta)

        logit Y D X1 X2 Z
        margins, dydx(D)
        global true_ame_mc = r(table)[1,1]
		
		*-----------------------------------------------------------------------
        * True rho_Y: partial correlation between Z and the outcome residual
        * from the SHORT model (omitting Z), after partialling Z on X1 X2.
        *
        * Step 1: fit the short logit (omitting Z) and compute Pearson residuals
        * Step 2: residualise Z on X1 X2 to get Zr
        * Step 3: standardise both and compute their correlation
        *
        * This is the population analogue of what the sensitivity command
        * computes in finite samples, so it gives the true rho_Y that the
        * benchmark analysis would recover if Z were observed.
        *-----------------------------------------------------------------------

        * Short model Pearson residuals
        logit Y D X1 X2
        predict double mu_short, pr
        gen double Yr_short = (Y - mu_short) / sqrt(mu_short*(1-mu_short))
        replace Yr_short = 0 if Yr_short == .
        sum Yr_short
        replace Yr_short = (Yr_short - r(mean)) / r(sd)

        * Treatment residual (standardised) — same as sensitivity command
        reg D X1 X2
        predict double Dr_short, resid
        sum Dr_short
        replace Dr_short = Dr_short / r(sd)

        * Residualise Z on X1 X2 (partialling out controls)
        reg Z X1 X2
        predict double Zr, resid
        sum Zr
        replace Zr = Zr / r(sd)

        * rho_D: partial correlation of Z with treatment residual
        corr Zr Dr_short
        global true_rhoD_mc = r(rho)

        * rho_Y: partial correlation of Z with outcome residual
        corr Zr Yr_short
        global true_rhoY_mc = r(rho)

        drop mu_short Yr_short Dr_short Zr
    }
    restore
end


capture program drop sim_v3
program define sim_v3, rclass

    syntax , nsim(integer) nobs(integer)

    *---------------------------------------------------------------------------
    * Fixed DGP parameters
    *---------------------------------------------------------------------------
    local g1   =  0.3
    local g2   =  0.2
    local intc = -0.5

    *---------------------------------------------------------------------------
    * Cell definitions
    *---------------------------------------------------------------------------
    local cell_beta  "0    0    0    0.5  0.5  0.5"
    local cell_delta "0.2  0.5  0.8  0.2  0.5  0.8"
    local cell_rho   "0.3  0.4  0.6  0.3  0.4  0.6"
    local cell_label `""null/very-weak" "null/moderate" "null/strong" "pos/very-weak" "pos/moderate" "pos/strong""'
    local ncells = 6

    *---------------------------------------------------------------------------
    * Results matrices: 6 cells x 3 methods
    *---------------------------------------------------------------------------
    tempname res_mae res_bias res_dir res_trueame
    matrix `res_mae'    = J(`ncells', 3, .)
    matrix `res_bias'   = J(`ncells', 3, .)
    matrix `res_dir'    = J(`ncells', 3, .)
    matrix `res_trueame'= J(`ncells', 1, .)
    matrix colnames `res_mae'  = CH_OLS CH_adapted Synthetic
    matrix colnames `res_bias' = CH_OLS CH_adapted Synthetic
    matrix colnames `res_dir'  = CH_OLS CH_adapted Synthetic
    matrix rownames `res_mae'  = c1 c2 c3 c4 c5 c6
    matrix rownames `res_bias' = c1 c2 c3 c4 c5 c6
    matrix rownames `res_dir'  = c1 c2 c3 c4 c5 c6

    *---------------------------------------------------------------------------
    * Cell loop
    *---------------------------------------------------------------------------
    forval cell = 1/`ncells' {

        local beta_d = real(word("`cell_beta'",  `cell'))
        local dlt    = real(word("`cell_delta'", `cell'))
        local rho0   = real(word("`cell_rho'",   `cell'))
        local lbl    : word `cell' of `cell_label'
        local null_cell = (`beta_d' == 0)

        di as text _n "{hline 70}"
        di as text "CELL `cell': `lbl'"
        di as text "  beta_D=`beta_d'  delta_Z=`dlt'  rho_D=`rho0'"
        di as text "{hline 70}"

        *-----------------------------------------------------------------------
        * True AME
        *-----------------------------------------------------------------------
        compute_true_ame, intc(`intc') beta(`beta_d') g1(`g1') g2(`g2') ///
            delta(`dlt') rho0(`rho0')
        local true_ame = $true_ame_mc
        matrix `res_trueame'[`cell',1] = `true_ame'
        di as text "  True AME = " %8.4f `true_ame'

        *-----------------------------------------------------------------------
        * Accumulators
        *-----------------------------------------------------------------------
        foreach m in ols cha syn {
            local `m'_sum_ae = 0
            local `m'_sum_se = 0
            local `m'_n_dir  = 0
            local `m'_valid  = 0
        }

        *-----------------------------------------------------------------------
        * Simulation loop
        *-----------------------------------------------------------------------
        forval sim = 1/`nsim' {
        quietly {

            *-------------------------------------------------------------------
            * 1. Generate data
            *-------------------------------------------------------------------
            clear
            set obs `nobs'
            gen X1 = rnormal()
            gen X2 = rnormal()
            gen D  = (0.3*X1 + rnormal() > 0)

            reg D X1 X2
            predict double Dr0, resid
            sum Dr0
            replace Dr0 = Dr0 / r(sd)
            gen double Z = `rho0'*Dr0 + sqrt(1-`rho0'^2)*rnormal()
            drop Dr0

            gen double eta = `intc' + `beta_d'*D + `g1'*X1 + `g2'*X2 + `dlt'*Z
            gen byte   Y   = runiform() < invlogit(eta)

            *-------------------------------------------------------------------
            * 2. Short logit
            *-------------------------------------------------------------------
            local skip = 0
            capture {
                logit Y D X1 X2
                margins, dydx(D)
                local me_short = r(table)[1,1]
                local se_short = r(table)[2,1]

                predict double mu_s, pr
                gen double Yr = (Y - mu_s) / sqrt(mu_s*(1-mu_s))
                replace Yr = 0 if Yr == .
                sum Yr
                replace Yr = (Yr - r(mean)) / r(sd)
                drop mu_s
            }
            if _rc != 0 local skip = 1
            if `skip' {
                drop _all
                set obs 1
                continue
            }

            * Treatment residual
            reg D X1 X2
            predict double Dr, resid
            sum Dr
            replace Dr = Dr / r(sd)

            corr Dr Yr
            local r_main = r(rho)

            *-------------------------------------------------------------------
            * 3. Benchmark partial correlations for Z
            *-------------------------------------------------------------------
            reg Z X1 X2
            predict double Zr, resid
            sum Zr
            replace Zr = Zr / r(sd)

            corr Zr Dr
            local rho_D_k = r(rho)
            corr Zr Yr
            local rho_Y_k = r(rho)

            *-------------------------------------------------------------------
            * METHOD 1: CH-OLS
            *-------------------------------------------------------------------
            local ols_adj = .
            capture {
                reg Y D X1 X2
                local b_ols  = _b[D]
                local se_ols = _se[D]
                local n_ols  = e(N)

                predict double Yr_ols, resid
                sum Yr_ols
                replace Yr_ols = (Yr_ols - r(mean)) / r(sd)

                reg D X1 X2
                predict double Dr_ols, resid
                sum Dr_ols
                replace Dr_ols = Dr_ols / r(sd)

                corr Zr Dr_ols
                local r2dz_ols = r(rho)^2
                corr Zr Yr_ols
                local r2yz_ols = r(rho)^2
                drop Yr_ols Dr_ols

                ch_ols_adjusted, b(`b_ols') se(`se_ols') n(`n_ols') ///
                    r2dz(`r2dz_ols') r2yz(`r2yz_ols') ncontrols(2)
                local ols_adj = $ch_b_adj
            }
            if _rc != 0 local ols_adj = .

            *-------------------------------------------------------------------
            * METHOD 2: CH-adapted
            *-------------------------------------------------------------------
            local cha_adj = .
            capture {
                local r2dz_cha = `rho_D_k'^2
                local r2yz_cha = `rho_Y_k'^2
                ch_ols_adjusted, b(`me_short') se(`se_short') n(`nobs') ///
                    r2dz(`r2dz_cha') r2yz(`r2yz_cha') ncontrols(2)
                local cha_adj = $ch_b_adj
            }
            if _rc != 0 local cha_adj = .

            *-------------------------------------------------------------------
            * METHOD 3: Synthetic confounder
            *-------------------------------------------------------------------
            local syn_adj = .
            local denom_syn = 1 - `r_main'^2
            if abs(`denom_syn') >= 1e-8 {
                local av  = (`rho_D_k' - `rho_Y_k'*`r_main') / `denom_syn'
                local bv  = (`rho_Y_k' - `rho_D_k'*`r_main') / `denom_syn'
                local csq = 1 - `av'^2 - `bv'^2 - 2*`av'*`bv'*`r_main'
                if `csq' >= -1e-10 {
                    local cv = sqrt(max(0, `csq'))
                    gen double eps_v = rnormal()
                    gen double Z_hyp = `av'*Dr + `bv'*Yr + `cv'*eps_v
                    capture {
                        logit Y D X1 X2 Z_hyp
                        margins, dydx(D)
                        local syn_adj = r(table)[1,1]
                    }
                    if _rc != 0 local syn_adj = .
                    drop Z_hyp eps_v
                }
            }

 *-------------------------------------------------------------------
            * 4. Accumulate errors
            * Direction error: adjustment moved AME in wrong direction
            * i.e. (adj - short) has opposite sign to (true - short)
            * Correct direction: adjustment should move AME toward true_ame
            *-------------------------------------------------------------------
            foreach m in ols cha syn {
                local adj = ``m'_adj'
                if `adj' != . {
                    local `m'_valid  = ``m'_valid'  + 1
                    local `m'_sum_ae = ``m'_sum_ae' + abs(`adj' - `true_ame')
                    local `m'_sum_se = ``m'_sum_se' + (`adj' - `true_ame')

                    * Direction of required adjustment
                    local required_dir = sign(`true_ame' - `me_short')
                    * Direction of actual adjustment
                    local actual_dir   = sign(`adj' - `me_short')

                    * Error if actual direction opposes required direction
                    * (ignore draws where required adjustment is exactly zero)
                    if `required_dir' != 0 {
                        if `actual_dir' != `required_dir' {
                            local `m'_n_dir = ``m'_n_dir' + 1
                        }
                    }
                }
            }

            drop Dr Yr Zr

        } // end quietly
        } // end sim loop

        *-----------------------------------------------------------------------
        * Store results
        *-----------------------------------------------------------------------
        local col = 1
        foreach m in ols cha syn {
            local dv = max(1, ``m'_valid')
            local mae_`m'  = ``m'_sum_ae' / `dv'
            local bias_`m' = ``m'_sum_se' / `dv'
            local dir_`m'  = 100 * ``m'_n_dir' / `dv'
            matrix `res_mae'[`cell',`col']  = `mae_`m''
            matrix `res_bias'[`cell',`col'] = `bias_`m''
            matrix `res_dir'[`cell',`col']  = `dir_`m''
            local col = `col' + 1
        }

        di as text "  Valid draws: ols=`ols_valid'  cha=`cha_valid'" ///
                   "  syn=`syn_valid'"
        di as text %32s " " %14s "CH-OLS" %14s "CH-adapted" %14s "Synthetic"
        di as text %32s "MAE" ///
            %14.4f `mae_ols'  %14.4f `mae_cha'  %14.4f `mae_syn'
        di as text %32s "Bias" ///
            %14.4f `bias_ols' %14.4f `bias_cha' %14.4f `bias_syn'
        if `null_cell' {
            di as text %32s "Made-worse rate (%)" ///
                %14.1f `dir_ols' %14.1f `dir_cha' %14.1f `dir_syn'
        }
        else {
            di as text %32s "Direction error rate (%)" ///
                %14.1f `dir_ols' %14.1f `dir_cha' %14.1f `dir_syn'
        }

    } // end cell loop

    *---------------------------------------------------------------------------
    * Summary tables
    *---------------------------------------------------------------------------

    * --- Table 1: MAE ---
    di as text _n "{hline 80}"
    di as text "TABLE 1: MEAN ABSOLUTE ERROR — lower is better"
    di as text "Prediction: MAE advantage of synthetic grows with delta_Z."
    di as text "{hline 80}"
    di as text %22s "Cell" %10s "delta_Z" %14s "CH-OLS" ///
               %14s "CH-adapted" %14s "Synthetic"
    di as text "{hline 80}"
    forval c = 1/`ncells' {
        local lbl   : word `c' of `cell_label'
        local dlt_c = real(word("`cell_delta'", `c'))
        di as result %22s "`lbl'" %10.1f `dlt_c' ///
            %14.4f `res_mae'[`c',1] ///
            %14.4f `res_mae'[`c',2] ///
            %14.4f `res_mae'[`c',3]
        * Separator between null and positive panels
        if `c' == 3 di as text %80s "{hline 80}"
    }
    di as text "{hline 80}"

    * --- Table 2: Bias ---
    di as text _n "{hline 80}"
    di as text "TABLE 2: BIAS (adjusted AME - true AME)"
    di as text "Positive=overcorrects, negative=undercorrects."
    di as text "Prediction: CH bias grows with delta_Z (O(delta^2))."
    di as text "{hline 80}"
    di as text %22s "Cell" %10s "delta_Z" %14s "CH-OLS" ///
               %14s "CH-adapted" %14s "Synthetic"
    di as text "{hline 80}"
    forval c = 1/`ncells' {
        local lbl   : word `c' of `cell_label'
        local dlt_c = real(word("`cell_delta'", `c'))
        di as result %22s "`lbl'" %10.1f `dlt_c' ///
            %14.4f `res_bias'[`c',1] ///
            %14.4f `res_bias'[`c',2] ///
            %14.4f `res_bias'[`c',3]
        if `c' == 3 di as text %80s "{hline 80}"
    }
    di as text "{hline 80}"

    * --- Table 3: Direction error ---
    di as text _n "{hline 80}"
    di as text "TABLE 3: DIRECTION ERROR RATE (%)"
di as text "% draws where adjustment moved AME away from true AME."
di as text "Applies to all cells consistently."
di as text "Should be near 0% for all methods when confounding is strong."
di as text "May be noisy when required adjustment is near zero (weak conf.)."
    di as text "{hline 80}"
    di as text %22s "Cell" %10s "delta_Z" %14s "CH-OLS" ///
               %14s "CH-adapted" %14s "Synthetic"
    di as text "{hline 80}"
    forval c = 1/`ncells' {
        local lbl   : word `c' of `cell_label'
        local dlt_c = real(word("`cell_delta'", `c'))
        di as result %22s "`lbl'" %10.1f `dlt_c' ///
            %14.1f `res_dir'[`c',1] ///
            %14.1f `res_dir'[`c',2] ///
            %14.1f `res_dir'[`c',3]
        if `c' == 3 di as text %80s "{hline 80}"
    }
    di as text "{hline 80}"

    * --- Ratio table: synthetic MAE advantage ---
    di as text _n "{hline 80}"
    di as text "TABLE 4: RATIO CH-OLS MAE / SYNTHETIC MAE"
    di as text "Values > 1 favour synthetic. Prediction: ratio grows with delta_Z."
    di as text "{hline 80}"
    di as text %22s "Cell" %10s "delta_Z" ///
               %24s "CH-OLS / Synthetic" %24s "CH-adapted / Synthetic"
    di as text "{hline 80}"
    forval c = 1/`ncells' {
        local lbl   : word `c' of `cell_label'
        local dlt_c = real(word("`cell_delta'", `c'))
        local ratio_ols = `res_mae'[`c',1] / max(1e-8, `res_mae'[`c',3])
        local ratio_cha = `res_mae'[`c',2] / max(1e-8, `res_mae'[`c',3])
        di as result %22s "`lbl'" %10.1f `dlt_c' ///
            %24.2f `ratio_ols' %24.2f `ratio_cha'
        if `c' == 3 di as text %80s "{hline 80}"
    }
    di as text "{hline 80}"
    di as text "Ratio = 1: methods equally accurate."
    di as text "Ratio > 1: synthetic more accurate by that factor."

    * --- True AMEs for reference ---
    di as text _n "{hline 80}"
    di as text "TRUE AMEs BY CELL (for reference)"
    di as text "{hline 80}"
    forval c = 1/`ncells' {
        local lbl : word `c' of `cell_label'
        di as result %22s "`lbl'" %10.4f `res_trueame'[`c',1]
    }
    di as text "{hline 80}"

    return matrix mae     = `res_mae'
    return matrix bias    = `res_bias'
    return matrix dir     = `res_dir'
    return matrix trueame = `res_trueame'

end


*===============================================================================
* Run
*===============================================================================
set seed 12345
sim_v3, nsim(10000) nobs(5000)


log close









*===============================================================================
* EMPIRICAL APPLICATION: Cognitive skills and labour supply at age 46
* BCS70 (1970 British Cohort Study)
*
* Treatment (D):   pc1 -- first principal component of five age-10
*                  cognitive test scores (Fun Maths Test; BAS word
*                  definitions, number recall, similarities, matrices;
*                  each age-standardised before the PCA). Constructed
*                  in datapreparation_2.do (line ~193).
* Outcomes (Y):    employedFT, unemployed, isManagerSupervisor -- all
*                  binary, from the 2016 wave (age 46), Study 8547.
* Controls:        perserverance_2_sd (teacher-rated perseverance,
*                  item j139) plus the full set1-set7 covariate blocks
*                  (birth, perinatal, age-5, age-10, medical, teacher's
*                  questionnaire, residual classroom-behaviour PCs).
* Benchmark:       c_mum_degree, c_dad_degree (parental degree-level
*                  qualification, age-5 parent questionnaire) -- the
*                  closest available analogue to the motheduc/fatheduc
*                  benchmarking used in the Mroz application, and
*                  already included as controls (set3), so the
*                  package's leave-one-out benchmark() logic applies
*                  to them directly without modification.
*
* IMPORTANT, FLAGGED EXPLICITLY RATHER THAN GLOSSED OVER:
*   1. This uses our own sens_* package (sensitivity_v17_validate_
*      full_with_sims.do), NOT senseweighter/senseweighter_hetero.
*   2. Restricted to the SINGLE 2016 (age 46) wave. The original
*      panel_labour_outcomes.dta is built via repeated 1:m merges
*      across six waves; our package has no clustering support, so
*      running it on the full stacked panel would treat repeated
*      observations on the same person as independent and understate
*      standard errors. Restricting to one wave avoids this rather
*      than papering over it.
*   3. hours and logwage/net_pay_week_95 are OUT OF SCOPE: our
*      `sensitivity` command only supports logit/probit/cloglog/
*      poisson/nbreg/ologit/mlogit -- no linear-model option exists,
*      so continuous labour-supply outcomes are not run here.
*   4. NO IPW/attrition weighting and NO cluster(bcsid): our package's
*      sens_fit calls plain logit with no [pweight=] support. This is
*      a real simplification relative to the existing IPW-weighted,
*      clustered pipeline in analysis_weights.do/cinelli_extended.do,
*      not an enhancement -- the results below should be read as an
*      illustration of the sensitivity method, not as a like-for-like
*      replacement for the team's existing, more careful estimates.
*===============================================================================

* ------------------------------------------------------------------------------
* Load the package (sens_abc, sens_build_z, sens_fit, sens_bisect, sensitivity, ...)
* ------------------------------------------------------------------------------
*do "/Users/user/Dropbox/Econometrics/sensitivity_v17_validate_full_with_sims.do"

* ------------------------------------------------------------------------------
* Load the covariates (run in this order: 3, 1, 2 -- per the team's own convention)
* ------------------------------------------------------------------------------

log using "/Users/user/Library/CloudStorage/Dropbox/Econometrics/politicalBehaviour/bcs70_sens.smcl", replace


qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"

keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)

* ------------------------------------------------------------------------------
* Redefine set1 to avoid the "i." operator, exactly as the team's existing
* sensitivity-section code already does (analysis_weights.do, line ~347):
* FWL-style residualisation inside our package's sens_abc/sensitivity machinery
* cannot consume i.varname the way some Stata commands can.
* ------------------------------------------------------------------------------
tab c_region_birth, gen(region_)
tab c_sexbirth, gen(sex_)
global set1 c_delivery_age_mum c_delivery_age_mum_2 region_1 region_2 region_3 region_4 region_5 region_6 region_7 region_8 region_9 region_10 c_mumsmoked c_abnormal_gest sex_2 c_birthweight

* ------------------------------------------------------------------------------
* Restrict to the single 2016 (age 46) wave to avoid the repeated-observations
* problem described above. panel_labour_outcomes.dta has one row per bcsid per
* wave; we keep only wave==2016 before merging, so the merge below is 1:1, not
* 1:m, and every row in the working dataset is one person.
* ------------------------------------------------------------------------------
preserve
    use "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta", clear
    keep if wave == 2016
    tempfile labour_age46
    save `labour_age46'
restore

merge 1:1 bcsid using `labour_age46'
drop if _merge == 2
drop _merge

* ------------------------------------------------------------------------------
* Full control list for this exercise.
* ------------------------------------------------------------------------------
global cogcontrols perserverance_2_sd $set1 $set2 $set3 $set4 $set5 $set6 $set7



* --- Profile plots for employedFT ---
sensitivity_profile employedFT pc1 $cogcontrols, model(logit) ///
    rhoy(0.1) ngrid(201) saving(profile_empFT_ry01) replace

sensitivity_profile employedFT pc1 $cogcontrols, model(logit) ///
    rhoy(0.3) ngrid(201) saving(profile_empFT_ry03) replace

* --- Profile plots for isManagerSupervisor ---
sensitivity_profile isManagerSupervisor pc1 $cogcontrols, model(logit) ///
    rhoy(0.1) ngrid(201) saving(profile_mgr_ry01) replace

sensitivity_profile isManagerSupervisor pc1 $cogcontrols, model(logit) ///
    rhoy(0.3) ngrid(201) saving(profile_mgr_ry03) replace
	
	
* ------------------------------------------------------------------------------
* Run the sensitivity package for each binary labour-supply outcome at age 46.
* benchmark(c_mum_degree c_dad_degree): both already sit in $set3 (controls),
* so the package's leave-one-out logic (drop the variable from controls,
* measure its own (rho_D, rho_Y), evaluate the implied bias) applies directly,
* mirroring how motheduc/fatheduc were benchmarked in the Mroz application.
* ------------------------------------------------------------------------------
local outcomes employedFT  isManagerSupervisor

foreach y of local outcomes {

    di as text _n "{hline 70}"
    di as text "COGNITIVE SKILLS (pc1) AND `y' AT AGE 46"
    di as text "{hline 70}"

    quietly count if !missing(`y', pc1)
    di as text "Non-missing observations for this outcome: " r(N)

    capture noisily sensitivity `y' pc1 $cogcontrols, model(logit)   ///
        rhod(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5)         ///
        rhoy(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5)         ///
        bisect signchange                                            ///
        benchmark( sex_2 c_motherStayhome10 c_fatherUnemployed10 pcB c_heightAt10_sd c_weightAt10_sd )

    if _rc != 0 {
        di as error "sensitivity failed for outcome `y' (rc=" _rc "). Likely"
        di as error "causes: separation/perfect prediction in this subsample,"
        di as error "or a control with zero variance after the wave restriction."
        di as error "Check `y' and the control set before re-running."
    }
}

logit employedFT  pc1 $cogcontrols
mfx
logit  isManagerSupervisor  pc1 $cogcontrols
mfx

log close
