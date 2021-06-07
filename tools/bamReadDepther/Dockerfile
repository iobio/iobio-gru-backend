FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl clang zlib1g-dev

COPY build.sh /
RUN chmod +x build.sh
RUN /build.sh
