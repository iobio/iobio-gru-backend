FROM ubuntu:18.04

RUN apt-get update && apt-get install -y curl git

ARG node_version=14.x
RUN curl -sL https://deb.nodesource.com/setup_${node_version} | bash - && \
    apt-get install -y nodejs

# We expect node binary to be named "node", but Ubuntu names it "nodejs", so we create a symbolic link
RUN ln -s /usr/bin/nodejs /usr/local/bin/node

RUN mkdir -p /iobio-gru-backend/tool_bin/lib/phenotypeExtractor

WORKDIR /iobio-gru-backend/tool_bin/lib/phenotypeExtractor

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3

COPY build.sh /
RUN /build.sh

COPY phenotypeExtractor.sh /
