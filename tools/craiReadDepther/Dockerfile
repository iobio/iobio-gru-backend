FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl git

COPY install_go.sh /
COPY build.sh /

RUN chmod +x install_go.sh
RUN /install_go.sh

ENV PATH "/usr/local/go/bin:${PATH}"

RUN chmod +x build.sh
RUN /build.sh
