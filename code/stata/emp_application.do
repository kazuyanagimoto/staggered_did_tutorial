************************************************
***   Staggered DID: Problems & Solutions    ***
***   Part II. Empirical Application         ***
***   Written by Yoshifumi Konishi           *** 
***   Date: 2023.3.15                        ***
************************************************     

/************ Empirical Settings ***************

Cicala (2022) "Imperfect Markets versus Imperfect Regulation in US Electricity
Generation" American Economic Review

- Effect of market dispatch on electricity generation costs & gains from trade
- DID design exploiting staggered rollouts of market dispatch
- Data: Hourly generation cost data at the generation-unit level 
    ==> Use collapsed daily data 

************************************************/

/***** Setting global *****/

clear
clear mata
clear matrix
set maxvar 10000
set matsize 10000

* Directory
global hdir = "E:\FLA075\Research\Staggered DID\日経学会2023\Cicala_AER_2022_Data"

/***** I. Overview of the Data *****/

cd "$hdir"
use "cicala_aer_2022_ready", clear

gen y = log_ideal_trade_surplus /* Gains from trade in log */
gen treat = market_operation /* Market dispatch */
gen tr_time = mkt_ym /* Treatment timing */
egen group = group(mkt_ym) /* Timing group */ 
gen time = ym(year,month) /* Time var */
egen pca_id = group(pca_abbrev99) /* PCA ID var */

* check distribution of treatment over time * 
gen temp = 1
bysort year month: egen obs_by_month = total(temp)
bysort year month: egen mkt_count = total(market_operation)
gen mkt_frac = mkt_count/obs_by_month
drop temp
bysort date: egen num_pca = count(pca_id)

* check how many timing groups? *
sum pca_id group
table treat, command(sum ideal_trade_surplus)

sort time
format time %tm
tw line mkt_frac time, ytitle("Fraction of Market Dispatch", margin(r=2)) ///
	xlabel(468 "99" 480 "00" 492 "01" 504 "02" 516 "03" ///
	528 "04" 540 "05" 552 "06" 564 "07" 576 "08" 588 "09" ///
	600 "10" 612 "11" 624 "12" 636 "13") xtitle(" ") 
graph save dispatch_timing, replace

/***** II. Diagnosis *****/

*** II-A. Goodman-Bacon Decomposition ***

xtset pca_id time
set matsize 10000
bacondecomp y treat year, ddetail gropt(msymbols(oh t))
* ddtiming y treat, i(id) t(time) savedata(nfd2) replace

* Memo: Goodman-Bacon does not allow for multiple obs. within a panel.
*       e.g., if time var (= treatment level) is monthly, observations cannot be
*       hourly or daily. Similarly, if id var (= treatment level) is state, 
*       observations cannot be at the household level.

twowayfeweights y pca_id time treat

*** II-B. Jakiela Diagnosis ***

*** prep weight & residualized outcomes ***
qui: reghdfe treat, a(year i.pca_id##i.month i.pca_id#c.log_load) res
predict tr_resid, resid
gen tr_resid2 = tr_resid^2
egen denom = total(tr_resid2)
gen w = tr_resid/denom

*** Diagnosis (1): Weight ***

tw (scatter group time if time < tr_time, ///
	msymbol(s) mlcolor(gs10) mfcolor(gs14)) ///
(scatter group time if time >= tr_time & tr_resid > 0, ///
	msymbol(s) mcolor(dkgreen*0.5)) ///
(scatter group time if time >= tr_time & tr_resid < 0, ///
	msymbol(s) mcolor(maroon)), ///
	aspect(0.5) plotregion(style(none)) ///
	ylabel(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6"7 "7" 8 "8" 9 "9" 10 "10" ///
	11 "11" 12 "12" 13 "13" 14 "14" 15 "15" 16 "16" 17 "17", ///
	angle(0) noticks labsize(small)) yscale(lstyle(none)) ///
	ytitle("Timing group", size(medsmall) margin(r=2)) ///
	legend(order(3 2 1) label(3 "Treatment observations - negative weight") ///
	label(2 "Treatment observations - positive weight") ///
	label(1 "Comparison observations") size(small) col(1)) ///
	xtitle(" ") xlabel(468 "99" 480 "00" 492 "01" 504 "02" 516 "03" ///
	528 "04" 540 "05" 552 "06" 564 "07" 576 "08" 588 "09" ///
	600 "10" 612 "11" 624 "12" 636 "13", labsize(small))
graph save jdiag_1, replace 

table group, stat(count treat) stat(mean treat) stat(sd treat) 
		
*** Diagnosis (2): Heterogeneity ***
qui: reghdfe y, a(year i.pca_id##i.month i.pca_id#c.log_load) res
predict y_resid, resid
tw (scatter y_resid tr_resid if treat == 0, m(oh) mc(maroon%40)) ///
	(scatter y_resid tr_resid if treat == 1, m(sh)  mc(dkgreen%40)) ///
	(lfit y_resid tr_resid if treat == 0, lc(maroon)) ///
	(lfit y_resid tr_resid if treat == 1, lc(dkgreen)) ///
	, xtitle("D residual") ytitle("Y residual") xlabel(-.5(.5).5) ///
	legend(order(1 "Treatment observations" 2 "Comparison observations")) ///
	legend(ring(0) pos(1) col(1) region(style(off))) 
graph save jdiag_2, replace  

* The following yield the same results *
reg y_resid tr_resid, cluster(pca_modate)
reghdfe y treat, absorb(datetime_region pca_month) cluster(pca_modate)

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
save "cicala_alt.dta", replace

*** 0. Benchmark TWFE ***
use "cicala_alt.dta", clear

reghdfe y lead* lag* [aweight=wgt], ///
	a(year i.pca_id##i.month i.pca_id#c.log_load) cluster(pca_modate)
* Note: lead1 coef is set to zero. 
*       Unlike in simulated data, all included leads/lags are estimated. 

event_plot, default_look stub_lag(lag#) stub_lead(lead#) together ///
	trimlead(12) trimlag(18) ///
    graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("Naive TWFE Event Study", size(medsmall) margin(b=3)) ///
	xlabel(-12(3)18) name(cg0, replace)) 
graph save cicala_g0, replace

*** 1. Fully-saturated TWFE ***

* (1-B) Sun and Abraham *
use "cicala_alt.dta", clear

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
graph save cicala_g1, replace

*** 2. Rolling methods ***

* (2-A) Callaway-Sant'Anna *
use "cicala_alt.dta", clear

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
graph save cicala_g2_2, replace

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
graph save cicala_g2, replace
*/

* (2-B) de Chaisemartin and D'Haultfoeuille *
use "cicala_alt.dta", clear

* Absorb FEs and Time-varying Covariates a la Caetano et al. (2022) *
qui reg y i.year i.pca_id##i.month i.pca_id#c.log_load ///
	if treat == 0, nocons /* weighting does not matter here */
predict ydc, residual

timer on 3
did_multiplegt ydc pca_id time treat, robust_dynamic dynamic(24) placebo(24) ///
	breps(50) seed(1) cluster(pca_modate) 
* Memo: Package allows for time-varying controls, but not factor variables.
*       Hence, residualize manually.	
timer off 3

event_plot e(estimates)#e(variances), default_look ///
	stub_lag(Effect_#) stub_lead(Placebo_#) together ///
	trimlead(12) trimlag(18) ///
	graph_opt(xtitle("Months since the event") ytitle("Estimates") ///
	title("de Chaisemartin and D'Haultfoeuille", size(medsmall) ///
	margin(b=3)) xlabel(-12(6)18) name(cg3, replace)) 
graph save cicala_g3, replace

*** 3. Imputation methods ***

* (3-A) Borusyak et al. *
use "cicala_alt.dta", clear

timer on 4
did_imputation y pca_id time tr_time [aweight=wgt], ///
	horizon(0/18) fe(year pca_id#month) unitc(log_load) ///
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
graph save cicala_g4_2, replace
	
* (3-B) Gardner *
use "cicala_alt.dta", clear

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
graph save cicala_g5, replace

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





 