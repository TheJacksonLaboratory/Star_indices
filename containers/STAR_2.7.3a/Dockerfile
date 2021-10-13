# continuumio/miniconda3:4.10.3-alpine
FROM continuumio/miniconda3@sha256:e007b47bbf7de9db9adec13588593d02cdaf0c32ee257f7a2f88737bf8b9945d
LABEL description="Dockerfile containing all the requirements for STAR RNA-seq read aligner" \
      author="ines@lifebit.ai"

ARG ENV_NAME="base"

RUN apk add --no-cache bash=5.1.4-r0 procps=3.3.17-r0

COPY environment.yml /
RUN conda env update -n ${ENV_NAME} -f environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

# Initialise bash for conda
RUN conda init bash

USER root
WORKDIR /data/
CMD ["bash"]
