FROM rocker/rstudio

RUN apt update && apt install -y \
    openssh-client libxt-dev\
    # Python
    python3 python3-pip

# R Packages
RUN R -e "install.packages(c('languageserver', 'renv'))"

# DVC Path
ENV PATH $PATH:~/.pip/bin