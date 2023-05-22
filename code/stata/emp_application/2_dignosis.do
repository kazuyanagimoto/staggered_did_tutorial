/***** II. Diagnosis *****/

use output/stata/emp_application/1_overview/cicala.dta, clear

*** II-A. Goodman-Bacon Decomposition ***

* xtset pca_id date
* bacondecomp y treat, ddetail gropt(msymbols(oh t))
* ddtiming y treat, i(id) t(time) savedata(nfd2) replace

* Memo: Goodman-Bacon does not allow for multiple obs. within a panel.
*       e.g., if time var (= treatment level) is monthly, observations cannot be
*       hourly or daily. Similarly, if id var (= treatment level) is state,
*       observations cannot be at the household level.

*** II-B. Jakiela Diagnosis ***

*** prep weight & residualized outcomes ***
qui: reghdfe treat, a(year i.pca_id##i.month i.pca_id#c.log_load) res
predict tr_resid if e(sample), resid
gen tr_resid2 = tr_resid^2
egen denom = total(tr_resid2)
gen w = tr_resid/denom
bysort group time: egen wm = mean(w)
* Note: Confirm distribution of weights over time is not uniform across groups
*       by running : scatter w time if group == `i' & time >= tr_time

*** Diagnosis (1): Weight ***

tw (scatter group time if time < tr_time, ///
	msymbol(s) mlcolor(gs10) mfcolor(gs14)) ///
(scatter group time if time >= tr_time & wm > 0, ///
	msymbol(s) mcolor(dkgreen*0.5)) ///
(scatter group time if time >= tr_time & wm < 0, ///
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
graph export output/stata/emp_application/2_dignosis/jdiag_1.png, replace

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
graph export output/stata/emp_application/2_dignosis/jdiag_2.png, replace

* The following yield the same results *
reg y_resid tr_resid, cluster(pca_modate)
reghdfe y treat, absorb(datetime_region pca_month) cluster(pca_modate)