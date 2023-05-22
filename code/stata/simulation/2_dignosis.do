/***** II. Diagnosis *****/

use output/stata/simulation/1_gen_data/data.dta, clear

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
bysort group time: egen wm = mean(w)
* Note: Confirm distribution of weights over time is not uniform across groups
*       by running : scatter w time if group == `i' & time >= tr_time

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

graph export output/stata/simulation/2_diagnosis/jakiela_weight.pdf, replace
graph export output/stata/simulation/2_diagnosis/jakiela_weight.png, replace

*** Diagnosis (2): Heterogeneity ***
qui: reg y i.time i.id
predict y_resid, resid

tw (scatter y_resid tr_resid if treat == 0, m(oh) mc(maroon%40)) ///
	(scatter y_resid tr_resid if treat == 1, m(sh)  mc(dkgreen%40)) ///
	(lfit y_resid tr_resid if treat == 0, lc(maroon)) ///
	(lfit y_resid tr_resid if treat == 1, lc(dkgreen)) ///
	, xtitle("D residual") ytitle("Y residual") ///
	title("Case 2. Heterogeneous/Dynamic Effect" ///
	"(Simulation 6 in Baker {it:et al.})", size(medsmall) margin(b=3))  ///
	legend(order(1 "Treatment observations" 2 "Comparison observations")) ///
	legend(ring(0) pos(1) col(1) region(style(off))) xsize(8.5)

graph export output/stata/simulation/2_diagnosis/jakiela_resid.pdf, replace
graph export output/stata/simulation/2_diagnosis/jakiela_resid.png, replace