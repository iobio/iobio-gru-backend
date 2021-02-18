#FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_102.0

USER root
WORKDIR /

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y \
        curl \
        git \
        xz-utils \
        python-pip \
        build-essential


ARG node_version=12.x
RUN curl -sL https://deb.nodesource.com/setup_${node_version} | bash - && \
    apt-get install -y nodejs

# We expect node binary to be named "node", but Ubuntu names it "nodejs", so we create a symbolic link
RUN ln -s /usr/bin/nodejs /usr/local/bin/node

ARG gru_version=dev
RUN git clone --depth 1 --branch ${gru_version} https://github.com/iobio/iobio-gru-backend

WORKDIR /iobio-gru-backend

RUN ["npm", "install"]

RUN ["mkdir", "tool_bin"]

# TODO: see if we can just copy what we need from the vep image rather than
# basing off it
#COPY --from=ensemblorg/ensembl-vep:latest /opt/vep /vep
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/samtools tool_bin/samtools-1.11
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/tabix tool_bin/tabix-1.11
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/bgzip tool_bin/bgzip-1.11
COPY --from=anderspitman/htslib:1.11 /bcftools-1.11/bcftools tool_bin/bcftools-1.11
COPY --from=anderspitman/bai-read-depther:latest /bamReadDepther tool_bin/baiReadDepther
COPY --from=anderspitman/crai-read-depther:latest /craiReadDepther tool_bin/craiReadDepther
COPY --from=anderspitman/bam-stats-alive:latest /bamstatsAlive/bamstatsAlive tool_bin/bamstatsAlive

RUN ln -s $PWD/tool_bin/bcftools-1.11 $PWD/tool_bin/bcftools
RUN ln -s $PWD/tool_bin/samtools-1.11 $PWD/tool_bin/samtools
RUN ln -s $PWD/tool_bin/tabix-1.11 $PWD/tool_bin/tabix
RUN ln -s $PWD/tool_bin/bgzip-1.11 $PWD/tool_bin/bgzip

EXPOSE 9001

ENTRYPOINT ["node", "src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin"]

