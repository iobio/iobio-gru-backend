FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl build-essential zlib1g-dev

COPY vcfReadDeptherHelper.sh /

COPY build.sh /
RUN /build.sh
