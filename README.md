# Staggered DiD Tutorial

## Tutorials

- [発表資料スライド](https://github.com/kazuyanagimoto/staggered_did_tutorial/blob/main/docs/slides/%E6%97%A5%E7%B5%8C%E5%AD%A6%E4%BC%9A%E3%83%81%E3%83%A5%E3%83%BC%E3%83%88%E3%83%AA%E3%82%A2%E3%83%AB_%E7%99%BA%E8%A1%A8%E8%B3%87%E6%96%99.pdf)
- [Stata演習スライド](https://github.com/kazuyanagimoto/staggered_did_tutorial/blob/main/docs/slides/%E6%97%A5%E7%B5%8C%E5%AD%A6%E4%BC%9A%E3%83%81%E3%83%A5%E3%83%BC%E3%83%88%E3%83%AA%E3%82%A2%E3%83%AB_Stata%E6%BC%94%E7%BF%92.pdf)
- [Comparison between R and Stata](https://github.com/kazuyanagimoto/staggered_did_tutorial/blob/main/docs/report/report.html)

## Replication (Stata & R)

### Stata

1. Clone this repository
1. Set the root directory in `setup.do` and run to install the packages
1. Run `simulation/1_gen_data.do`-`simulation/3_estimation.do`
1. Run `emp_application/1_overviw.do`-`emp_application/3_estimation.do`

*Important Caveat*

In the empirical part, we experienced that the computation of `csdid` (Callaway and Sant'Anna) crushed
the computer system (not only the Stata application) with a RAM of 32GB or less. 
We strongly recommend you run it only if your computer has 64GB or larger RAM.

### R

1. Clone this repository
1. Open the local repository as an R project
    (`staggered_did_tutorial.Rproj`)
1. Run `renv::restore()` in the R console to install packages
1. Run `simulation/1_gen_data.R`-`simulation/3_estimation.R`
1. Run `emp_application/1_overview.R`-`emp_application/3_estimation.R`


For docker-users, between the step 1 and 2, implement either

- **Rstudio:** `docker-compose up --build` and open `localhost:8787` in Browser
- **VSCode:** Dev Containers "Rebuild and Reopen in Container"



## License and Credits

Copyright (c) 2023 Yoshifumi Konishi, Kazuharu Yanagimoto, and Hayato Umetani.

Tutorial slides and Stata codes are written by Konishi, and R codes are written by Yanagimoto and Umetani.
Source codes and data are under the Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) License, see [LICENSE](https://github.com/kazuyanagimoto/staggered_did_tutorial/blob/main/LICENSE).
Note that the original hourly-level data of `data/cicala_aer_2022_ready.dta` are provided by
Cicala (2022) under the CC BY-NC 4.0 License.
The data is converted to a daily-level and cleaned by Konishi.

## References

Cicala, Steve. Hourly U.S. Electricity Load. Nashville, TN: American Economic Association [publisher], 2022. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2022-01-29. https://doi.org/10.3886/E146801V1