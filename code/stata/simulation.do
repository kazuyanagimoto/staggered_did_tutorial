************************************************
***   Staggered DID: Problems & Solutions    ***
***   Part I. Basics using Simulated Data    ***
***   Written by Yoshifumi Konishi           *** 
***   Date: 2023.3.14                        ***
************************************************     


/***** Install Packages *****/

/*
*** (1) Utility packages ***
ssc install reghdfe, replace
ssc install ftools, replace
ssc install event_plot, replace
ssc install addplot, replace
ssc install esttab, replace

*** (2) Decomposition packages ***
ssc install bacondecomp, replace /* Goodman-Bacon decomp */
ssc install ddtiming, replace /* (Old) Goodman-Bacon decomp */
ssc install eventstudyweights, replace /* Sun-Abraham decomp */
ssc install twowayfeweights, replace /* deCDH decomp */

*** (3) DID packages ***
ssc install stackedev, replace /* Stacked regression */ 
ssc install eventstudyinteract, replace /* Sun-Abraham*/
scc install csdid, replace /* Callaway-Sant'Anna */
ssc install did_multiplegt, replace /* deCDH */
ssc install did_imputation, replace /* Borusyak et al. */
ssc install did2s, replace /* Gardner */

*/

/***** I. Generating Data *****/

*** Case 2. Heterogeneous, Dynamic Effect ***
clear 
set obs 1000
gen id = _n
set seed 100000
gen ai = rnormal(0,0.5)
forvalues t = 0(1)35 {
scalar t = rnormal(0,0.5)
gen lt`t' = t
}
gen rn = runiform()
gen group = 1 if rn <= 17/50
replace group = 2 if rn > 17/50
replace group = 3 if rn > 35/50

gen tr_time = 1989 if group == 1
replace tr_time = 1998 if group == 2
replace tr_time = 2007 if group == 3 

reshape long lt, i(id) j(time)
replace time = time + 1980
gen d89 = (time >= 1989)
gen d98 = (time >= 1998)
gen d07 = (time >= 2007) 

gen treat = d89 if group == 1
replace treat = d98 if group == 2
replace treat = d07 if group ==3 

gen y = ai + lt + rnormal(0,0.5)
replace y = y + rnormal(0.5,0.2)*(time-1989)*treat if group == 1
replace y = y + rnormal(0.3,0.2)*(time-1998)*treat if group == 2
replace y = y + rnormal(0.1,0.2)*(time-2007)*treat if group == 3

gen window = 1 if time >= 1989 & time <= 1997
replace window = 2 if time >= 1997 & time <= 2007
replace window = 3 if time >= 2007 & time <= 2015

gen effect =.
replace effect = (0.5)*(time-1989) if group == 1 & time >= 1989
replace effect = (0.3)*(time-1998) if group == 2 & time >= 1998
replace effect = (0.1)*(time-2007) if group == 3 & time >= 2007
bysort group window: sum effect
sum effect

bysort group time: egen ymean = mean(y)
bysort group time: egen ymax = max(y)
bysort group time: egen ymin = min(y)
replace ymax = 20 if ymax > 20

twoway (rarea ymax ymin time if group == 1, lc(black%0) fc(gray%5)) ///
(rarea ymax ymin time if group == 2, lc(black%0) fc(gray%5)) ///
(rarea ymax ymin time if group == 3, lc(black%0) fc(gray%5)) ///
(line ymean time if group == 1, lw(medthick) lc(dkgreen)) ///
(line ymean time if group == 2, lw(medthick) lc(navy)) ///
(line ymean time if group == 3, lw(medthick) lc(maroon)) ///
, xsize(5) legend(order(4 "1989" 5 "1998" 6 "2007") ///
region(lp(blank)) symxs(*0.5) symys(*0.5) size(medsmall) ///
ring(0) pos(11) rows(3)) xlabel(1980(10)2015) ///
xline(1989.5 1998.5 2007.5, lp(dash) lc(black) lw(thin)) ///
title("Case 2. Heterogeneous/Dynamic Effect""(Simulation 6 in Baker {it:et al.})", ///
size(medsmall) margin(b=3)) yscale(r(0 20)) ylabel(-5(5)20) xtitle(" ") name(g2, replace) 

* TWFE regression *
reghdfe y treat, absorb(id time) cl(id)
* xtset id time
* xtreg y treat i.time, fe vce(cl id)

/***** II. Diagnosis *****/

*** II-A. Goodman-Bacon Decomposition ***

xtset id time
bacondecomp y treat, msymbols(oh t) /* Current version, aggregates early vs. late */
* ddtiming y treat, i(id) t(time) savedata(nfd2) replace /* Old version, each timing comparison */
* do "sim_bacondecomp.do"

* Memo: Similar packages by Sun-Abraham and deCDH, but Jakiela diagnosis
*       is more intuitive and easier to customize for quick diagnosis. 

*** II-B. Jakiela Diagnosis ***

*** prep weight & residualized outcomes ***
qui: reg treat i.time i.id
predict tr_resid, resid
gen tr_resid2 = tr_resid^2
egen denom = total(tr_resid2)
gen w = tr_resid/denom

*** Diagnosis (1): Weight ***
tw (scatter group time if time < tr_time, ms(s) mlc(gs10) mfc(gs14) msize(vlarge)) ///
(scatter group time if time >= tr_time & tr_resid > 0, ms(s) mc(dkgreen*0.5) msize(vlarge)) ///
(scatter group time if time >= tr_time & tr_resid < 0, ms(s) mc(maroon) msize(vlarge)), ///
	aspect(0.1) plotregion(style(none)) ///
	ylabel(1 "Treat in 1989" 2 "Treat in 1998" 3 "Treat in 2007" , ///
		angle(0) noticks labsize(small)) ytitle(" " ) yscale(lstyle(none)) ///
		legend(order(3 2 1) label(3 "Treatment observations - negative weight") ///
		label(2 "Treatment observations - positive weight") ///
		label(1 "Comparison observations") size(small) col(1)) xtitle(" ") xlabel(1980(5)2015)
		
*** Diagnosis (2): Heterogeneity ***
qui: reg y i.time i.id
predict y_resid, resid

tw (scatter y_resid tr_resid if treat == 0, m(oh) mc(maroon%40)) ///
(scatter y_resid tr_resid if treat == 1, m(sh)  mc(dkgreen%40)) ///
(lfit y_resid tr_resid if treat == 0, lc(maroon)) ///
(lfit y_resid tr_resid if treat == 1, lc(dkgreen)) ///
, xtitle("D residual") ytitle("Y residual") ///
legend(order(1 "Treatment observations" 2 "Comparison observations")) ///
legend(ring(0) pos(1) col(1) region(style(off))) xsize(8.5)

/***** III. Alternative Estimators *****/

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
save "sim2_full.dta", replace

*** 0. Benchmark TWFE ***
use "sim2_full.dta", clear

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
use "sim2_full.dta", clear

keep if (time >= 1989 - 5) & (time <= 1989 + 5)
save "sim2_sub1.dta", replace

use "sim2_full.dta", clear
keep if (time >= 1998 - 5) & (time <= 1998 + 5)
keep if group >= 2
save "sim2_sub2.dta", replace

use "sim2_full.dta", clear
keep if (time >= 2007 - 5) & (time <= 2007 + 5)
keep if group == 3

append using "sim2_sub1.dta"
append using "sim2_sub2.dta"

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

use "sim2_full.dta", clear

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
use "sim2_full.dta", clear

forvalues g = 2/3 {
	forvalues l = 1/27 {
	gen lead`l'_g`g' = lead`l'*(group == `g')
	}
}
forvalues g = 2/3 {
	forvalues l = 0/27 {
	gen lag`l'_g`g' = lag`l'*(group == `g')
	}
}

reghdfe y lag* lead*, a(time group) cluster(id)

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
	reghdfe y lag* lead*, a(time group) cluster(id)
	margins, expression(_b[lag`i']+_b[lag`i'_g2]*(2/3)+_b[lag`i'_g3]*(1/3)) post
	mat blag`i' = e(b)
	mat vlag`i' = e(V)
	replace beta = blag`i'[1,1] if grid == `i'+6
	replace var = vlag`i'[1,1] if grid == `i'+6
}
forvalues i = 1/5 {
	reghdfe y lag* lead*, a(time group) cluster(id)
	margins, expression(_b[lead`i']+_b[lead`i'_g2]*(2/3)+_b[lead`i'_g3]*(1/3)) post
	mat blead`i' = e(b)
	mat vlead`i' = e(V)
	replace beta = blead`i'[1,1] if grid == 6-`i'
	replace var = vlead`i'[1,1] if grid == 6-`i'
}

gen lower = beta - 0.96*sqrt(var) in 1/11
gen upper = beta + 0.96*sqrt(var) in 1/11
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
use "sim2_full.dta", clear

timer on 2
csdid y, ivar(id) time(time) gvar(tr_time) notyet ///
    agg(event) wboot driwp cluster(id) 
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
use "sim2_full.dta", clear

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
use "sim2_full.dta", clear

timer on 4
did_imputation y id time tr_time, autosample horizon(0/5) pretrend(5)
timer off 4

event_plot, default_look stub_lag(tau#) stub_lead(pre#) together trimlead(5) trimlag(5) ///
    graph_opt(xtitle("Periods since the event") ytitle("Estimates") ///
	title("Borusyak, Jaravel, Spiess", size(medsmall) margin(b=3))  ///
	xlabel(-5(1)5) name(g6, replace)) 
addplot: (scatteri 0 0 2 5, xlabel(-5(1)5) recast(line) lp(dash) lc(maroon))
	
* (3-B) Gardner *
use "sim2_full.dta", clear

timer on 5
reg y i.id i.time if treat == 0, nocons
predict yhat, residual
reg yhat lead* lag* if d89 != 1 | time < 2007, nocons cluster(id) /// drop negative weight group
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
graph combine g4 g5 g6 g7, rows(2) xsize(5.5)
