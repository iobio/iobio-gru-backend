FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl git xz-utils bzip2 libfontconfig1

COPY build.sh /
RUN chmod +x build.sh
RUN /build.sh

COPY clinReport.sh /clin-report
