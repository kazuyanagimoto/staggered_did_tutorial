************************************************
***   Staggered DID: Problems & Solutions    ***
***   Part II. Empirical Application         ***
***   Written by Yoshifumi Konishi           *** 
***   Date: 2023.5.22                        ***
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

/***** I. Overview of the Data *****/

use "data/cicala_aer_2022_ready.dta", clear

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
graph export output/stata/emp_application/1_overview/dispatch_timing.pdf, replace
graph export output/stata/emp_application/1_overview/dispatch_timing.png, replace

save output/stata/emp_application/1_overview/cicala.dta, replace