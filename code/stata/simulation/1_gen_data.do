global RootDir "path-to-local-directry"
cd $RootDir

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

graph export output/stata/simulation/1_gen_data/baker22_sim6.pdf, replace

* TWFE regression *
reghdfe y treat, absorb(id time) cl(id)
* xtset id time
* xtreg y treat i.time, fe vce(cl id)

save output/stata/simulation/1_gen_data/data.dta, replace
