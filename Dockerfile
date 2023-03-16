FROM rocker/rstudio

RUN apt update && apt install -y \
    openssh-client libxt-dev tree

# R Packages
RUN R -e "install.packages(c('languageserver', 'renv'))"
