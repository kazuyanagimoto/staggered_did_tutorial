*** (1) Utility packages ***
ssc install reghdfe, replace
ssc install avar, replace /* eventsutydyinteract dependency */
ssc install ftools, replace
ssc install event_plot, replace
ssc install addplot, replace
ssc install estout, replace

*** (2) Decomposition packages ***
ssc install bacondecomp, replace /* Goodman-Bacon decomp */
ssc install eventstudyweights, replace /* Sun-Abraham decomp */
ssc install twowayfeweights, replace /* deCDH decomp */

*** (3) DID packages ***
ssc install stackedev, replace /* Stacked regression */ 
ssc install eventstudyinteract, replace /* Sun-Abraham*/
ssc install csdid, replace /* Callaway-Sant'Anna */
ssc install did_multiplegt, replace /* deCDH */
ssc install did_imputation, replace /* Borusyak et al. */
ssc install did2s, replace /* Gardner */


* Directory

global projDir = "PATH/TO/YOUR/DIRECTORY"

if "`c(username)'" == "kazuharu" {
  global projDir = "\\wsl.localhost\Ubuntu\home\kazuharu\github\staggered_did_tutorial"  
}
if "`c(username)'" == "umetanihayato" {
  global projDir = "/Users/umetanihayato/git/staggered_did_tutorial"
}

cd "$projDir"