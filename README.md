# Staggered DiD Tutorial

本Githubサイトは、2023年5月28日に日本経済学会春季大会にて開催されたチュートリアル・セッション（共催：日本学術会議 数量的経済・政策分析分科会）「DIDの計量経済手法の近年の展開」のサポートサイトとしてスライド資料（講義編・演習基礎編・演習応用編）、Stataコード、Rコードを提供しています．

## Tutorials

- [発表資料スライド](https://kazuyanagimoto.com/staggered_did_tutorial/slides/SDID_Tutorial_Main.pdf)
- [Stata演習スライド](https://kazuyanagimoto.com/staggered_did_tutorial/slides/SDID_Tutorial_Stata.pdf)
- [Comparison between R and Stata](https://kazuyanagimoto.com/staggered_did_tutorial/code/stata_vs_r.html)

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
1. By default, benchmarks are not calculated. Set `run_benchmark <- TRUE` for it.

For docker-users, between the step 1 and 2, implement either

- **Rstudio:** `docker-compose up --build` and open `localhost:8787` in Browser
- **VSCode:** Dev Containers "Rebuild and Reopen in Container"

## License and Credits

スライド資料およびStataコードは小西祥文（慶應義塾大学）が作成，Rコードは研究補助として柳本和春さん(CEMFI)と梅谷隼人さん(神戸大学)に作成頂きました．教育・研究目的での利用は引用の上, 商用利用の場合は私（小西祥文）の許諾を得てからご利用下さい．

Tutorial slides and Stata codes are written by Yoshifumi Konishi (Keio Univ.), and R codes are written by Kazuharu Yanagimoto (CEMFI) and Hayato Umetani (Kobe Univ.). Source codes and data are under the Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) License, see [LICENSE](https://github.com/kazuyanagimoto/staggered_did_tutorial/blob/main/LICENSE).

Note: The data set used in Stata/R demonstrations (`cicala_aer_2022_ready.dta`) comes from Cicala (AER, 2022) under the CC BY-NC 4.0 License. The original data are hourly observations but are converted to daily-level observations by Yoshifumi Konishi.

## References

Cicala, Steve. Hourly U.S. Electricity Load. Nashville, TN: American Economic Association [publisher], 2022. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2022-01-29. https://doi.org/10.3886/E146801V1
