# This is intended to run in Local Development (dev) and Github Actions (test/prod)
# BUILD_ENV options (dev, test, prod) dev for local testing and test for github actions testing on prod ready code
ARG BUILD_ENV="prod"
ARG MAINTAINER="kimn@ssi.dk;"
ARG BIFROST_COMPONENT_NAME="bifrost_chewbbaca"
ARG FORCE_DOWNLOAD=true

#---------------------------------------------------------------------------------------------------
# Programs for all environments
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.12.0 as build_base
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD ARG BUILD_ENV
ONBUILD ARG MAINTAINER
ONBUILD LABEL \
    BIFROST_COMPONENT_NAME=${BIFROST_COMPONENT_NAME} \
    description="Docker environment for ${BIFROST_COMPONENT_NAME}" \
    environment="${BUILD_ENV}" \
    maintainer="${MAINTAINER}"
ONBUILD RUN \
    conda install -yq -c conda-forge -c bioconda -c default snakemake-minimal==5.7.1;


#---------------------------------------------------------------------------------------------------
# Base for dev environement
#---------------------------------------------------------------------------------------------------
FROM build_base as build_dev
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD COPY /components/${BIFROST_COMPONENT_NAME} /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY /lib/bifrostlib /bifrost/lib/bifrostlib
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}/
ONBUILD RUN \
    pip install -r requirements.txt; \
    pip install --no-cache -e file:///bifrost/lib/bifrostlib; \
    pip install --no-cache -e file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for production environment
#---------------------------------------------------------------------------------------------------
FROM build_base as build_prod
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY ./ ./
ONBUILD RUN \
    pip install -e file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for test environment (prod with tests)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
FROM build_base as build_test
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY ./ ./
ONBUILD RUN \
    pip install -r requirements.txt \
    pip install -e file:///bifrost/components/${BIFROST_COMPONENT_NAME}/



#---------------------------------------------------------------------------------------------------
# Additional resources
#---------------------------------------------------------------------------------------------------
FROM build_${BUILD_ENV}
ARG BIFROST_COMPONENT_NAME
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
RUN \
    apt-get update &&  apt-get install -y -q --fix-missing \
        build-essential \
        zlib1g-dev \
        libmagic-dev \
        nano \
        less; \
    conda install -c bioconda blast=2.9; \
    conda install -c bioconda prodigal=2.6.3; \
    git clone https://github.com/B-UMMI/chewBBACA.git /opt/chewbbaca; \
    pip install -r /opt/chewbbaca/CHEWBBACA/requirements.txt; \
    pip install -e file:///opt/chewbbaca;
    #chewBBACA.py DownloadSchema -sp 4 -sc 1 -o resources/Campylobacter_jejuni_ns;
    #chewBBACA.py DownloadSchema -sp 5 -sc 1 -o resources/Escherichia_coli_ns; \
    #chewBBACA.py DownloadSchema -sp 6 -sc 1 -o resources/Listeria_monocytogenes_ns; \
    #chewBBACA.py DownloadSchema -sp 8 -sc 1 -o resources/Salmonella_enterica_ns;

#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENTRYPOINT ["python3", "-m", "bifrost_chewbbaca"]
CMD ["python3", "-m", "bifrost_chewbbaca", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------
