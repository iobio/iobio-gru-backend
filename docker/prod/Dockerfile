FROM ubuntu:18.04

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

COPY --from=anderspitman/htslib:1.11 /samtools-1.11/samtools tool_bin/samtools-1.11
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/tabix tool_bin/tabix-1.11
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/bgzip tool_bin/bgzip-1.11
COPY --from=anderspitman/bai-read-depther:latest /bamReadDepther tool_bin/baiReadDepther
COPY --from=anderspitman/bam-stats-alive:latest /bamstatsAlive/bamstatsAlive tool_bin/bamstatsAlive

EXPOSE 9001

ENTRYPOINT ["node", "src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin"]
