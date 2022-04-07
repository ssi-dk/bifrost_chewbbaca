# This is intended to run in Local Development (dev) and Github Actions (test/prod)
# BUILD_ENV options (dev, test, prod) dev for local testing and test for github actions testing on prod ready code
ARG BUILD_ENV="prod"
ARG MAINTAINER="kimn@ssi.dk;"
ARG BIFROST_COMPONENT_NAME="bifrost_chewbbaca"
ARG FORCE_DOWNLOAD=true

#---------------------------------------------------------------------------------------------------
# Programs for all environments
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.10.3 as build_base
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
    # For 'make' needed for kma
    apt-get update &&  apt-get install -y -q --fix-missing \
        build-essential \
        zlib1g-dev \
        libmagic-dev \
        nano \
        less; \
    conda install -c bioconda blast=2.12.0; \
    conda install -c bioconda prodigal=2.6.3; \
    conda install -c anaconda beautifulsoup4; \
    pip install -q \
        python-dateutil==2.8.1 \
        chewbbaca==2.8.5; \
    python html_download.py -a https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/ -o resources/Escherichia_cgMLSTv1; \
    python html_download.py -a https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/ -o resources/Salmonella_cgMLSTv2; \
    python html_download.py -a https://enterobase.warwick.ac.uk/schemes/Yersinia.cgMLSTv1/ -o resources/Yersinia_cgMLSTv1; \
    python download_alleles_bigsdb.py --base_url http://rest.pubmlst.org --database pubmlst_campylobacter_seqdef --scheme_id 4 --dir resources/Campylobacter_cgMLST; \
    python download_alleles_bigsdb.py --base_url https://bigsdb.pasteur.fr/api --database pubmlst_t21listeria_seqdef --scheme_id 3 --dir resources/Listeria_cgMLST; \
    for dir in resources/*/; do if ls $dir/*.gz 1> /dev/null 2>&1; then gunzip $dir/*.gz; fi; done; \
    chewBBACA.py PrepExternalSchema -i resources/Listeria_cgMLST/ -o resources/Listeria_cgMLST_prepped/ --ptf /opt/conda/lib/python3.9/site-packages/CHEWBBACA/prodigal_training_files/Listeria_monocytogenes.trn --cpu 4; \
    rm -r resources/Listeria_cgMLST; \
    chewBBACA.py PrepExternalSchema -i resources/Escherichia_cgMLSTv1/ -o resources/Escherichia_cgMLSTv1_prepped/ --ptf /opt/conda/lib/python3.9/site-packages/CHEWBBACA/prodigal_training_files/Escherichia_coli.trn --cpu 4; \
    rm -r resources/Escherichia_cgMLSTv1;
#    conda install python=3.6.5; \
#    conda install -c bioconda chewbbaca=2.0.16

#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENTRYPOINT ["python3", "-m", "bifrost_chewbbaca"]
CMD ["python3", "-m", "bifrost_chewbbaca", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------
