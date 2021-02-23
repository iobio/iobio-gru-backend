FROM ubuntu:18.04

RUN apt-get update && apt-get install -y git curl build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libgsl-dev

COPY build.sh /
RUN chmod +x build.sh
RUN /build.sh

COPY gtEnricher.sh /
RUN chmod +x /gtEnricher.sh
