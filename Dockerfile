FROM rocker/rstudio

RUN apt update && apt install -y --no-install-recommends\
    openssh-client libxt-dev

# R Packages
RUN R -e "install.packages(c('languageserver', 'renv'))"
