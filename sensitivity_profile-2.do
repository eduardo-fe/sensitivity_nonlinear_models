*===============================================================================
* sensitivity_profile
*
* Standalone diagnostic tool for the synthetic confounder sensitivity
* analysis. For a FIXED rho_Y, sweeps rho_D across the feasible interval
* on a fine grid, computes the adjusted AME and t-statistic at each point,
* and produces a two-panel graph showing:
*   - Top panel: adjusted AME as a function of rho_D, with a horizontal
*     line at AME = 0 (sign-change reference).
*   - Bottom panel: |t|-statistic as a function of rho_D, with a
*     horizontal line at z_{alpha/2} (explain-away reference).
*
* Threshold crossings are visible as the points where the AME curve
* crosses zero (sign change) and where the |t| curve crosses the
* critical value (explain away).
*
* REQUIRES: sens_abc, sens_build_z, sens_fit from the sensitivity package.
*   Run/load the main sensitivity_v18.do file first, or paste the three
*   helper programs (sens_abc, sens_build_z, sens_fit) before this file.
*
* SYNTAX:
*   sensitivity_profile Y D [X1 X2 ...], model(string) rhoy(real)
*       [ngrid(integer 201) seed(integer 12345) level(real 95)
*        baseoutcome(integer 1) outcome(string) atmeans
*        saving(string) replace scanmargin(real 0.005)]
*
* EXAMPLES:
*   sensitivity_profile Y D X1 X2, model(logit) rhoy(0.3)
*   sensitivity_profile Y D X1 X2, model(logit) rhoy(0.3) ngrid(501)
*   sensitivity_profile Y D X1 X2, model(mlogit) rhoy(0.3) ///
*       baseoutcome(1) outcome(2)
*   sensitivity_profile Y D X1 X2, model(logit) rhoy(0.3) ///
*       saving(profile_ry03) replace
*
*===============================================================================

capture program drop sensitivity_profile
program define sensitivity_profile, rclass
    syntax varlist(min=2) [if] [in],  ///
        Model(string)                 ///
        RHOY(real)                    ///
        [ NGRid(integer 201)          ///
          seed(integer 12345)         ///
          ATMeans                     ///
          BASEoutcome(integer 1)      ///
          Level(real 95)              ///
          outcome(string)             ///
          saving(string)              ///
          replace                     ///
          SCANMargin(real 0.005) ]

    marksample touse

    gettoken Y varlist : varlist
    gettoken D rest    : varlist
    local X "`rest'"

    local model = lower("`model'")
    if !inlist("`model'","logit","probit","cloglog","poisson", ///
                         "nbreg","ologit","mlogit") {
        di as error "Unsupported model: `model'"
        exit 198
    }

    if "`model'" == "mlogit" & "`outcome'" == "" {
        di as error "mlogit requires outcome() option"
        exit 198
    }

    global sens_X "`X'"

    local alpha = 1 - `level'/100
    local z     = invnormal(1 - `alpha'/2)
    local atm   = 0
    if "`atmeans'" != "" local atm = 1

    set seed `seed'

    di as text _n "{hline 70}"
    di as text "SENSITIVITY PROFILE: `model'"
    di as text "{hline 70}"
    di as text "Treatment: `D'"
    di as text "Outcome:   `Y'"
    if "`X'" != "" di as text "Controls:  `X'"
    di as text "Fixed rho_Y = " %6.3f `rhoy'
    di as text "Grid points: `ngrid'"
    di as text "Confidence level: `level'%"
    if "`outcome'" != "" di as text "MNL outcome: `outcome'"

    *===========================================================================
    * STEP 1: Outcome residuals from short model
    *===========================================================================

    if "`model'" == "logit" {
        qui logit `Y' `D' `X' if `touse'
        tempvar mu Yr
        qui predict double `mu' if `touse', pr
        qui gen double `Yr' = (`Y' - `mu') / sqrt(`mu'*(1-`mu')) if `touse'
        qui replace `Yr' = 0 if `Yr' == . & `touse'
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
        qui replace `Yr' = 0 if `Yr' == . & `touse'
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
        tempvar pr_o Yr
        qui predict double `pr_o' if `touse', outcome(`outcome')
        qui gen double `Yr' = (`Y' == real("`outcome'")) - `pr_o' if `touse'

        * Short-model AME for this outcome
        if `atm' == 1 qui margins, dydx(`D') atmeans predict(outcome(`outcome'))
        else          qui margins, dydx(`D') predict(outcome(`outcome'))
        local me_short = r(table)[1,1]
        local se_short = r(table)[2,1]
    }

    * Standardise outcome residual
    qui sum `Yr' if `touse'
    if r(sd) < 1e-8 {
        di as error "Outcome residual SD near zero"
        exit 498
    }
    qui replace `Yr' = (`Yr' - r(mean)) / r(sd) if `touse'

    *===========================================================================
    * STEP 2: Treatment residual (standardised)
    *===========================================================================
    tempvar Dr
    qui reg `D' `X' if `touse'
    qui predict double `Dr' if `touse', resid
    qui sum `Dr' if `touse'
    if r(sd) < 1e-8 {
        di as error "Treatment residual SD near zero"
        exit 498
    }
    qui replace `Dr' = `Dr' / r(sd) if `touse'

    *===========================================================================
    * STEP 3: Short-model AME and residual correlation
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
    }

    qui corr `Dr' `Yr' if `touse'
    local r = r(rho)

    di as text _n "Short model AME: " %9.4f `me_short' ///
        " (SE: " %7.4f `se_short' ")"
    di as text "r = cor(D_resid, Y_resid) = " %7.4f `r'
    di as text "t_short = " %7.2f (`me_short'/`se_short')

    *===========================================================================
    * STEP 4: Feasible interval
    *===========================================================================
    local disc_feas = (1 - (`r')^2) * (1 - (`rhoy')^2)
    if `disc_feas' < 0 {
        di as error "No feasible region for rho_Y = `rhoy'"
        exit 498
    }

    local rd_lo = max(-0.999, `r'*`rhoy' - sqrt(`disc_feas'))
    local rd_hi = min( 0.999, `r'*`rhoy' + sqrt(`disc_feas'))

    * Nudge if boundary is infeasible
    sens_abc, rd(`rd_hi') ry(`rhoy') rho(`r')
    if $sens_feasible == 0 local rd_hi = `rd_hi' - 0.001
    sens_abc, rd(`rd_lo') ry(`rhoy') rho(`r')
    if $sens_feasible == 0 local rd_lo = `rd_lo' + 0.001

    * Apply scanmargin
    local feas_width = `rd_hi' - `rd_lo'
    local rd_lo = `rd_lo' + `scanmargin' * `feas_width'
    local rd_hi = `rd_hi' - `scanmargin' * `feas_width'

    di as text "Feasible rho_D range: [" %6.3f `rd_lo' ", " %6.3f `rd_hi' "]"

    *===========================================================================
    * STEP 5: Fixed epsilon
    *===========================================================================
    tempvar eps_fixed
    qui gen double `eps_fixed' = rnormal() if `touse'
    global sens_eps "`eps_fixed'"

    *===========================================================================
    * STEP 6: Fine-grid evaluation
    *===========================================================================
    di as text _n "Evaluating `ngrid' grid points..."

    * Store results in a temporary dataset
    tempname results
    matrix `results' = J(`ngrid', 5, .)
    * col 1 = rho_D, col 2 = ME, col 3 = SE, col 4 = t-stat, col 5 = success

    local n_success = 0
    local n_fail    = 0

    forval g = 1/`ngrid' {
        local rd_g = `rd_lo' + (`g'-1) * (`rd_hi' - `rd_lo') / (`ngrid' - 1)

        matrix `results'[`g', 1] = `rd_g'

        sens_abc, rd(`rd_g') ry(`rhoy') rho(`r')
        if $sens_feasible == 0 {
            matrix `results'[`g', 5] = 0
            local n_fail = `n_fail' + 1
            continue
        }

        sens_build_z, dresid(`Dr') yresid(`Yr') touse(`touse')
        sens_fit, model(`model') yvar(`Y') dvar(`D') touse(`touse') ///
                  atm(`atm') baseoutcome(`baseoutcome') outcome("`outcome'")
        capture drop _Z_hyp_

        if $sens_failed == 0 & $sens_me != . {
            matrix `results'[`g', 2] = $sens_me
            matrix `results'[`g', 3] = $sens_se
            matrix `results'[`g', 4] = $sens_me / $sens_se
            matrix `results'[`g', 5] = 1
            local n_success = `n_success' + 1
        }
        else {
            matrix `results'[`g', 5] = 0
            local n_fail = `n_fail' + 1
        }

        * Progress indicator every 50 points
        if mod(`g', 50) == 0 {
            di as text "  ... `g' / `ngrid' done"
        }
    }

    di as text "  Done. Successful: `n_success' / `ngrid'"
    if `n_fail' > 0 {
        di as text "  Failed evaluations: `n_fail'"
    }

    *===========================================================================
    * STEP 7: Build plot dataset and graph
    *===========================================================================

    * Preserve current data, build plot dataset
    preserve
    clear
    qui set obs `ngrid'

    qui gen double rho_D  = .
    qui gen double ME     = .
    qui gen double SE     = .
    qui gen double tstat  = .
    qui gen double abs_t  = .
    qui gen byte   ok     = 0

    forval g = 1/`ngrid' {
        qui replace rho_D = `results'[`g', 1] in `g'
        if `results'[`g', 5] == 1 {
            qui replace ME    = `results'[`g', 2] in `g'
            qui replace SE    = `results'[`g', 3] in `g'
            qui replace tstat = `results'[`g', 4] in `g'
            qui replace abs_t = abs(`results'[`g', 4]) in `g'
            qui replace ok    = 1                      in `g'
        }
    }

    * Drop failed points for clean plotting
    qui drop if ok == 0

    * Confidence band for AME
    qui gen double ME_lo = ME - `z' * SE
    qui gen double ME_hi = ME + `z' * SE

    * Formatting
    local ry_fmt : di %5.2f `rhoy'
    local me_fmt : di %6.3f `me_short'
    local t_fmt  : di %5.1f (`me_short'/`se_short')
    local out_label ""
    if "`outcome'" != "" local out_label ", outcome `outcome'"

    * ---- Panel 1: AME ----
    twoway ///
        (rarea ME_lo ME_hi rho_D, color(gs14) lwidth(none)) ///
        (line ME rho_D, lcolor(navy) lwidth(medthick)) ///
        (function y = 0, range(rho_D) lcolor(cranberry) lpattern(dash) lwidth(medium)) ///
        (function y = `me_short', range(rho_D) lcolor(gs8) lpattern(shortdash) lwidth(thin)) ///
        , ///
        title("Adjusted AME as a function of {&rho}{sub:D}" ///
              "({&rho}{sub:Y} = `ry_fmt'`out_label')", size(medium)) ///
        xtitle("{&rho}{sub:D}") ///
        ytitle("Adjusted AME") ///
        ylabel(, angle(horizontal) format(%7.3f)) ///
        xlabel(, format(%5.2f)) ///
      	legend(order(2 "Adjusted AME" 1 "`level'% CI" ///
       3 "AME = 0" 4 "Short-model AME") ///
       rows(1) size(small) position(6) lstyle(none) region(lcolor(white))) ///
        scheme(s2color) ///
		graphregion(color(white)) plotregion(color(white)) ///
        name(sens_profile_ame, replace)

    * ---- Panel 2: |t|-statistic ----
    twoway ///
        (line abs_t rho_D, lcolor(navy) lwidth(medthick)) ///
        (function y = `z', range(rho_D) lcolor(cranberry) lpattern(dash) lwidth(medium)) ///
        , ///
        title("|t|-statistic as a function of {&rho}{sub:D}" ///
              "({&rho}{sub:Y} = `ry_fmt'`out_label')", size(medium)) ///
        xtitle("{&rho}{sub:D}") ///
        ytitle("|t|") ///
        ylabel(, angle(horizontal) format(%5.1f)) ///
        xlabel(, format(%5.2f)) ///
        yline(`z', lcolor(cranberry) lpattern(dash) lwidth(medium)) ///
        legend(order(1 "|t|-statistic" 2 "z = `z'") ///
               rows(1) size(small) position(6) lstyle(none) region(lcolor(white)) ) ///
         graphregion(color(white)) plotregion(color(white)) ///
        scheme(s2mono) ///
        name(sens_profile_tstat, replace)

    * ---- Combined ----
    graph combine sens_profile_ame sens_profile_tstat, ///
        cols(1) iscale(0.8) ///
        title("Sensitivity Profile: `model'", size(medium)) ///
        subtitle("Short-model AME = `me_fmt' (t = `t_fmt')" ///
                 "Fixed {&rho}{sub:Y} = `ry_fmt', grid = `ngrid' points", ///
                 size(small)) ///
        name(sens_profile_combined, replace)

    * ---- Save if requested ----
    if "`saving'" != "" {
        graph export "`saving'_combined.png", name(sens_profile_combined) ///
            width(2400) `replace'
        graph export "`saving'_ame.png", name(sens_profile_ame) ///
            width(2400) `replace'
        graph export "`saving'_tstat.png", name(sens_profile_tstat) ///
            width(2400) `replace'
        di as text _n "Graphs saved to `saving'_combined.png, " ///
            "`saving'_ame.png, `saving'_tstat.png"
    }

    * ---- Save dataset if requested ----
    if "`saving'" != "" {
        qui save "`saving'_data.dta", `replace'
        di as text "Data saved to `saving'_data.dta"
    }

    *===========================================================================
    * STEP 8: Locate and report crossings
    *===========================================================================
    di as text _n "{hline 70}"
    di as text "THRESHOLD CROSSINGS DETECTED"
    di as text "{hline 70}"

    local sign0 = sign(`me_short')

    * Explain-away crossings: where |t| crosses z
    local n_ea = 0
    local prev_above = .
    local prev_rd    = .
    local N = _N
    forval i = 1/`N' {
        local this_t    = abs_t[`i']
        local this_rd   = rho_D[`i']
        local this_above = (`this_t' >= `z')

        if `prev_above' != . & `this_above' != `prev_above' {
            local n_ea = `n_ea' + 1
            local cross_rd = (`prev_rd' + `this_rd') / 2
            di as text "  Explain-away crossing #`n_ea' at rho_D ~ " ///
                %7.4f `cross_rd'
        }
        local prev_above = `this_above'
        local prev_rd    = `this_rd'
    }
    if `n_ea' == 0 {
        di as text "  No explain-away crossings (effect remains significant" ///
            " throughout)"
    }

    * Sign-change crossings: where ME crosses zero
    local n_sc = 0
    local prev_sign = .
    local prev_rd   = .
    forval i = 1/`N' {
        local this_me   = ME[`i']
        local this_rd   = rho_D[`i']
        local this_sign = sign(`this_me')

        if `prev_sign' != . & `this_sign' != `prev_sign' & `this_sign' != 0 {
            local n_sc = `n_sc' + 1
            local cross_rd = (`prev_rd' + `this_rd') / 2
            di as text "  Sign-change crossing #`n_sc' at rho_D ~ " ///
                %7.4f `cross_rd'
        }
        local prev_sign = `this_sign'
        local prev_rd   = `this_rd'
    }
    if `n_sc' == 0 {
        di as text "  No sign-change crossings (sign never reverses)"
    }
    di as text "{hline 70}"

    restore

    * Return values
    return scalar me_short = `me_short'
    return scalar se_short = `se_short'
    return scalar r        = `r'
    return scalar rho_y    = `rhoy'
    return scalar rd_lo    = `rd_lo'
    return scalar rd_hi    = `rd_hi'
    return scalar n_ea     = `n_ea'
    return scalar n_sc     = `n_sc'
    return matrix profile  = `results'

    * Clean up
    macro drop sens_X sens_eps sens_a sens_b sens_c sens_feasible ///
               sens_me sens_se sens_failed sens_delta

end




* Generate the example data
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
gen eta2  = -0.5 + 0.5*D + 0.2*X1 + 0.3*X2 + 0.4*Z
gen eta3  = -1.0 + 0.7*D + 0.1*X1 + 0.4*X2 + 0.5*Z
gen denom = 1 + exp(eta2) + exp(eta3)
gen u     = runiform()
gen Y_mul = cond(u < 1/denom, 1, cond(u < (1+exp(eta2))/denom, 2, 3))

* Now run the profiles
sensitivity_profile Y_bin D X1 X2, model(logit) rhoy(0.3)
sensitivity_profile Y_mul D X1 X2, model(mlogit) rhoy(0.5) baseoutcome(1) outcome(1) saving("/Users/user/Library/CloudStorage/Dropbox/Econometrics/politicalBehaviour/outcome2_rho05") replace
