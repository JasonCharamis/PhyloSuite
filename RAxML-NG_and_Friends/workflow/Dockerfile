## Dockerfile to build Docker image for this Snakemake workflow.

FROM condaforge/mambaforge:latest
LABEL maintainer = "Jason Charamis" 
LABEL contact    = "jason.charamis@gmail.com"
LABEL build_date = "2023-11-30"
LABEL version    = "v.0.0.1-dev"

# Step 1: Install the essentials

RUN apt-get update &&\
    apt-get install -f -y &&\
    apt-get install build-essential -y &&\
    apt-get install git -y

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/phylo.yaml
#   prefix: /conda-envs/e7731bb49f40c7aec2d59fb33378cd57
#   name: phylo
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - mafft
#       - trimal
#       - pypythia
#       - modeltest-ng
#       - raxml-ng

RUN mkdir -p /conda-envs/e7731bb49f40c7aec2d59fb33378cd57
COPY envs/phylo.yaml /conda-envs/e7731bb49f40c7aec2d59fb33378cd57/environment.yaml

# Step 3: Generate conda environments

RUN mamba env create --prefix /conda-envs/e7731bb49f40c7aec2d59fb33378cd57 --file /conda-envs/e7731bb49f40c7aec2d59fb33378cd57/environment.yaml && \
    mamba clean --all -y

# Step 4: Install and test Snakemake

RUN /conda-envs/e7731bb49f40c7aec2d59fb33378cd57/bin/pip install wheel
RUN /conda-envs/e7731bb49f40c7aec2d59fb33378cd57/bin/pip install snakemake

ENV PATH="/conda-envs/e7731bb49f40c7aec2d59fb33378cd57/bin:$PATH"

# Step 5: Download repository and put scripts on PATH
RUN git clone https://github.com/JasonCharamis/Snakemake-workflow-for-RAxML-NG-based-phylogenetic-analysis.git

ENV PATH="Snakemake-workflow-for-RAxML-NG-based-phylogenetic-analysis/workflow/scripts:$PATH"

CMD ["Everything installed!"]
