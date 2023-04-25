/***** III. Alternative Estimators *****/

use output/stata/simulation/1_gen_data/data.dta, clear

gen rel_time = time - tr_time
sum tr_time
gen lastcohort = (tr_time == r(max)) // dummy for the latest- or never-treated cohort
forvalues l = 0/27 {
	gen lag`l' = (rel_time == `l')
}
forvalues l = 1/27 {
	gen lead`l' = (rel_time == -`l')
}
replace lead1 = 0 // normalize first_treat = -1 to be zero
save "output/stata/simulation/3_estimation/sim2_full.dta", replace

*** 0. Benchmark TWFE ***
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

reghdfe y lead* lag* , a(id time) cluster(id)
* Note 1: lead1 coef is set to zero. This will drop lags 18-27.
* Note 2: Runnin reghdfe y lag* lead* instead will drop leads 19-27 and lag 27,
*         and it will get us unbiased estimates since it "happens" to drop wrong comparisons.

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together ///
    graph_opt(xtitle("Periods since the event") ytitle("Estimates") ///
	title("Naive TWFE Event Study", size(medsmall) margin(b=3)) xlabel(-5(1)5) ///
	name(g0, replace)) trimlead(5) trimlag(5)
addplot: (scatteri 0 -5 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))


*** 1. Fully-saturated TWFE ***

* (1-A) Stacked regression *

* subsetting & stacking *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

keep if (time >= 1989 - 5) & (time <= 1989 + 5)
save "output/stata/simulation/3_estimation/sim2_sub1.dta", replace

use "output/stata/simulation/3_estimation/sim2_full.dta", clear
keep if (time >= 1998 - 5) & (time <= 1998 + 5)
keep if group >= 2
save "output/stata/simulation/3_estimation/sim2_sub2.dta", replace

use "output/stata/simulation/3_estimation/sim2_full.dta", clear
keep if (time >= 2007 - 5) & (time <= 2007 + 5)
keep if group == 3

append using "output/stata/simulation/3_estimation/sim2_sub1.dta"
append using "output/stata/simulation/3_estimation/sim2_sub2.dta"

reghdfe y lag* lead* , a(id time) cluster(id)

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together trimlead(5) trimlag(5) ///
    graph_opt(xtitle("Periods since the event") ytitle("Estimates") xlabel(-5(1)5) ///
	title("Stacked Regression", size(medsmall) margin(b=3)) xlabel(-5(1)5) ///
	name(g1, replace))
	addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

* Note 1: "stackedev" package requires never-treated units.
*         Hence cannot be used in this example b/c all units eventually receive treatment.
* Note 2: Runnding reghdfe y lead* lag* removes "correct" comparions due to multicoliniarity.

* (1-B) Sun and Abraham *

use "output/stata/simulation/3_estimation/sim2_full.dta", clear

timer clear
timer on 1
eventstudyinteract y lead* lag*, vce(cluster i) absorb(id time) cohort(tr_time) ///
	control_cohort(lastcohort)
timer off 1

event_plot e(b_iw)#e(V_iw), default_look stub_lag(lag#) stub_lead(lead#) ///
	together trimlead(5) trimlag(5) ///
	graph_opt(xtitle("Periods since the event") ytitle("Estimates") xlabel(-5(1)5) ///
	title("Sun and Abraham", size(medsmall) margin(b=3)) name(g2, replace))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

/*
eventstudyweights y lead* lag*, vce(cluster i) absorb(id time) ///
	cohort(tr_time) rel_time(rel_time) control_cohort(lastcohort) saveweights("weights")
import excel "weights.xlsx", clear firstrow
keep L1event first_treat rel_time
reshape wide L1event, i(rel_time) j(first_treat)
twoway line L1event* rel_time, xline(1, lp(dash) lw(thin) lc(black)) ///
   yline(0, lp(dash) lw(thin) lc(black)) xlabel(-30(5)30) ///
   ytitle("weights") xtitle("time relative to treatment")  ///
   title("Sun-Abraham Decomp""(Weights for t = 1)", size(medsmall) margin(b=2))
*/

* (1-C) Wooldridge *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

forvalues g = 1/2 {
	forvalues l = 1/27 {
	gen lead`l'_g`g' = lead`l'*(group == `g')
	}
}
forvalues g = 1/2 {
	forvalues l = 0/27 {
	gen lag`l'_g`g' = lag`l'*(group == `g')
	}
}

reghdfe y lag*_* lead*_*, a(time group) cluster(id)

/*
lincom (lag0 + lag0_g2*(2/3) + lag0_g3*(1/3))
lincom (lag1 + lag1_g2*(2/3) + lag1_g3*(1/3))
lincom (lag2 + lag2_g2*(2/3) + lag2_g3*(1/3))
lincom (lag3 + lag3_g2*(2/3) + lag3_g3*(1/3))
lincom (lag4 + lag4_g2*(2/3) + lag4_g3*(1/3))
lincom (lag5 + lag5_g2*(2/3) + lag5_g3*(1/3))
*/

gen grid = _n in 1/11
gen beta = .
gen var = .

forvalues i = 0/5 {
	reghdfe y lag*_* lead*_*, a(time group) cluster(id)
	margins, expression((_b[lag`i'_g1]+_b[lag`i'_g2])/2) post
	mat blag`i' = e(b)
	mat vlag`i' = e(V)
	replace beta = blag`i'[1,1] if grid == `i'+6
	replace var = vlag`i'[1,1] if grid == `i'+6
}
forvalues i = 1/5 {
	reghdfe y lag*_* lead*_*, a(time group) cluster(id)
	margins, expression(_b[lead`i'_g1]+_b[lead`i'_g2]/2) post
	mat blead`i' = e(b)
	mat vlead`i' = e(V)
	replace beta = blead`i'[1,1] if grid == 6-`i'
	replace var = vlead`i'[1,1] if grid == 6-`i'
}

gen lower = beta - 1.96*sqrt(var) in 1/11
gen upper = beta + 1.96*sqrt(var) in 1/11
replace grid = grid - 6
tw (rarea upper lower grid, lc(navy) fc(navy) fi(30)) ///
(connected beta grid, lc(navy) mc(navy)) ///
, xtitle("Periods since the event") ytitle("Estimates")  ///
	title("Wooldridge", size(medsmall) margin(b=3)) xlabel(-5(1)5) ///
	name(g3, replace) legend(off) ///
	yline(0, lc(black*0.5)) xline(0, lc(black*0.5) lp(dash))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

*** 2. Rolling methods ***

* (2-A) Callaway-Sant'Anna *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

timer on 2
csdid y, ivar(id) time(time) gvar(tr_time) notyet ///
    agg(event) wboot driwp
	* Memo: Default uses "uniform" SEs to account for multiple testing.
timer off 2
estat event, estore(cs)

event_plot cs, default_look stub_lag(Tp#) stub_lead(Tm#) together ///
	trimlead(5) trimlag(5) ///
	graph_opt(xtitle("Periods since the event") ytitle("Estimates")  ///
	title("Callaway and Sant'Anna", size(medsmall) margin(b=3)) xlabel(-5(1)5) ///
	name(g4, replace))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

* (2-B) de Chaisemartin and D'Haultfoeuille *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

timer on 3
did_multiplegt y id time treat, robust_dynamic ///
	dynamic(5) placebo(5) breps(100) cluster(id)
timer off 3

event_plot e(estimates)#e(variances), default_look ///
	stub_lag(Effect_#) stub_lead(Placebo_#) together ///
	trimlead(5) trimlag(5) ///
	graph_opt(xtitle("Periods since the event") ytitle("Estimates") ///
	title("de Chaisemartin and D'Haultfoeuille", size(medsmall) margin(b=3)) ///
	xlabel(-5(1)5) name(g5, replace))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

/*
twowayfeweights y id time treat, type(feTR) path("dCDH_decomp.dta")
use "dCDH_decomp.dta", clear
tw (scatter weight Time_TWFE)
*/

*** 3. Imputation methods ***

* (3-A) Borusyak et al. *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

timer on 4
did_imputation y id time tr_time, autosample horizon(0/5) pretrend(5)
timer off 4

event_plot, default_look stub_lag(tau#) stub_lead(pre#) together trimlead(5) trimlag(5) ///
    graph_opt(xtitle("Periods since the event") ytitle("Estimates") ///
	title("Borusyak, Jaravel, Spiess", size(medsmall) margin(b=3))  ///
	xlabel(-5(1)5) name(g6, replace))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

* (3-B) Gardner *
use "output/stata/simulation/3_estimation/sim2_full.dta", clear

timer on 5
reg y i.id i.time if treat == 0, nocons
predict yhat, residual
reg yhat lead* lag* if d89 != 1 | time < 2007, nocons cluster(id) // drop negative weight group
timer off 5

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together trimlead(5) trimlag(5) ///
    graph_opt(xtitle("Periods since the event") ytitle("Estimates") xlabel(-5(1)5) ///
	title("Gardner", size(medsmall) margin(b=3)) xlabel(-5(1)5) ///
	name(g7, replace))
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))

* Note : The following Stata package yeilds the same result.
*        did2s y, first_stage(i.id i.time) second_stage(lead* lag*) ///
*              treatment(treat) cluster(id)

*** Combine All Results ****
forvalues i = 1/5 {
	timer list `i'
}
graph combine g0 g1 g2 g3, rows(2) xsize(5.5)
graph export output/stata/simulation/3_estimation/alt_est1.pdf, replace
graph export output/stata/simulation/3_estimation/alt_est1.png, replace
graph combine g4 g5 g6 g7, rows(2) xsize(5.5)
graph export output/stata/simulation/3_estimation/alt_est2.pdf, replace
graph export output/stata/simulation/3_estimation/alt_est2.png, replace
