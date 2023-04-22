clear
clear mata
clear matrix
set maxvar 10000
use output/stata/emp_application/1_overview/cicala.dta, clear

/***** III. Alternative Estimators *****/

* Define relative time *

* Memo: Cicala (2022) estimate only wihtin 24-month window.
gen rel_time = time - tr_time
sum tr_time
gen lastcohort = (tr_time == r(max)) // dummy for the latest-treated cohort
gen never = mi(mkt_ym) // dummy for never-treated cohort
forvalues l = 0/24 {
	gen lag`l' = (rel_time == `l')
} 
forvalues l = 1/24 {
	gen lead`l' = (rel_time == -`l')
}
replace lead1 = 0 // normalize first_treat = -1 to be zero	
save "output/stata/emp_application/3_estimation/cicala_alt.dta", replace

*** 0. Benchmark TWFE ***
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

reghdfe y lead* lag* [aweight=wgt], ///
	a(year i.pca_id##i.month i.pca_id#c.log_load) cluster(pca_modate)
* Note: lead1 coef is set to zero. 
*       Unlike in simulated data, all included leads/lags are estimated. 

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together ///
	trimlead(12) trimlag(18) ///
    graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("Naive TWFE Event Study", size(medsmall) margin(b=3)) ///
	xlabel(-12(3)18) name(cg0, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g0, replace

*** 1. Fully-saturated TWFE ***

* (1-B) Sun and Abraham *
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

timer clear
timer on 1
eventstudyinteract y lead* lag* [aweight=wgt], vce(cluster pca_modate) ///
	absorb(year i.pca_id##i.month i.pca_id#c.log_load) cohort(tr_time) ///
	control_cohort(never)
timer off 1

event_plot e(b_iw)#e(V_iw), default_look stub_lag(lag#) stub_lead(lead#) ///
	together trimlead(12) trimlag(18) ///
	graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("Sun and Abraham", size(medsmall) margin(b=3)) xlabel(-12(3)18) ///
	name(cg1, replace))
graph save output/stata/emp_application/3_estimation/cicala_g1, replace

*** 2. Rolling methods ***

* (2-A) Callaway-Sant'Anna *
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

* Make sure to code "never treated group" = 0 *
replace tr_time = 0 if mi(tr_time) 

* Absorb FEs and Time-varying Covariates a la Caetano et al. (2022) *
qui reg y i.year i.pca_id##i.month i.pca_id#c.log_load ///
	if treat == 0, nocons /* weighting does not matter here */
predict ycs, residual

timer on 2
csdid ycs, time(time) gvar(tr_time) agg(event) ///
    wboot(reps(50)) rseed(1) cluster(pca_modate) 
* Memo 1: Drop ivar() to allow for multiple obs within a panel.
* Memo 2: Use wgt as covariate to weight obs.
timer off 2
* estat event, estore(cs)

event_plot, default_look stub_lag(Tp#) stub_lead(Tm#) together ///
	trimlead(12) trimlag(18) ///
	graph_opt(xtitle("Months since the event") ytitle("Estimates")  ///
	title("Callaway and Sant'Anna", size(medsmall) margin(b=3)) ///
	xlabel(-12(6)18) name(cg2, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g2, replace

/*
timer on 2
csdid ycs, time(time) gvar(tr_time) agg(event) ///
    wboot(reps(50)) rseed(1) long replace 
* Memo 1: Drop ivar() to allow for multiple obs within a panel.
* Memo 2: Use wgt as covariate to weight obs.
timer off 2
estat event, estore(cs)

event_plot cs, default_look stub_lag(Tp#) stub_lead(Tm#) together ///
	trimlead(12) trimlag(18) ///
	graph_opt(xtitle("Months since the event") ytitle("Estimates")  ///
	title("Callaway and Sant'Anna", size(medsmall) margin(b=3)) ///
	xlabel(-12(6)18) name(cg2, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g2, replace
*/

* (2-B) de Chaisemartin and D'Haultfoeuille *
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

* Absorb FEs and Time-varying Covariates a la Caetano et al. (2022) *
qui reg y i.year i.pca_id##i.month i.pca_id#c.log_load ///
	if treat == 0, nocons /* weighting does not matter here */
predict ydc, residual

timer on 3
did_multiplegt ydc pca_id time treat, robust_dynamic dynamic(18) placebo(12) ///
	breps(50) seed(1) cluster(pca_modate) firstdiff_placebo weight(wgt)
* Memo: Package allows for time-varying controls, but not factor variables.
*       Hence, residualize manually.	
timer off 3

event_plot e(estimates)#e(variances), default_look ///
	stub_lag(Effect_#) stub_lead(Placebo_#) together ///
	trimlead(12) trimlag(18) ///
	graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("de Chaisemartin and D'Haultfoeuille", size(medsmall) ///
	margin(b=3)) xlabel(-12(6)18) name(cg3, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g3, replace

*** 3. Imputation methods ***

* (3-A) Borusyak et al. *
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

timer on 4
did_imputation y pca_id time tr_time [aweight=wgt], ///
	horizon(0/18) fe(year month pca_id#month) unitc(log_load) ///
	autosample pretrend(12)
	* Full model get the error message: Could not run imputation 
	* for some observations b/c some absorbed variables/FEs 
	* are collinear in the D==0 subsample but in the full sample.
timer off 4

event_plot, default_look stub_lag(tau#) stub_lead(pre#) together ///
	trimlead(12) trimlag(18) ///
    graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("Borusyak, Jaravel, Spiess", size(medsmall) margin(b=3)) ///
	xlabel(-12(6)18) ///
	name(cg4, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g4, replace

	
* (3-B) Gardner *
use "output/stata/emp_application/3_estimation/cicala_alt.dta", clear

timer on 5
qui reg y i.year i.pca_id##i.month i.pca_id#c.log_load [aweight=wgt] if treat == 0, nocons
predict yhat, residual
reg yhat lead* lag* [aweight=wgt], nocons cluster(pca_modate)
timer off 5

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together ///
	trimlead(12) trimlag(18) ///
    graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("Gardner", size(medsmall) margin(b=3)) xlabel(-12(6)18) ///
	name(cg5, replace)) 
graph save output/stata/emp_application/3_estimation/cicala_g5, replace


* Note : The following Stata package yeilds the same result.
*        did2s y, first_stage(i.id i.time) second_stage(lead* lag*) treatment(treat) cluster(id)

*** Combine All Results ****
forvalues i = 1/5 {
	timer list `i'
}
forvalues i = 0/5 {
	graph use cicala_g`i', name(cg`i',replace) 
}
graph combine cg0 cg1 cg2 cg3, rows(3) xsize(5.5)
graph combine cg0 cg1 cg4 cg5, rows(3) xsize(5.5)


