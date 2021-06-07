FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl git cmake build-essential zlib1g-dev

COPY build.sh /
RUN chmod +x build.sh
RUN /build.sh
