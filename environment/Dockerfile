# hash:sha256:7261d71e0df0cf9bf521413efdaa5299fde887a8c434ab5660ce78cac3868d1a
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.12.4-R4.3.3-JupyterLab4.2.4-RStudiorstudio-server-2024.04.2-764-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN pip install -U --no-cache-dir \
    anndata \
    h5py \
    numpy \
    pandas \
    scipy

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("shiny")' \
    && Rscript -e 'remotes::install_version("bslib")' \
    && Rscript -e 'remotes::install_version("arrow")' \
    && Rscript -e 'remotes::install_version("dplyr")' \
    && Rscript -e 'remotes::install_version("ggplot2")' \
    && Rscript -e 'remotes::install_version("scales")' \
    && Rscript -e 'remotes::install_version("shinycssloaders")' \
    && Rscript -e 'remotes::install_version("cachem")' \
    && Rscript -e 'remotes::install_version("DT")'
