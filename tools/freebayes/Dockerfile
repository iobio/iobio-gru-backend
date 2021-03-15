FROM ubuntu:18.04

RUN apt-get update && apt-get -y install curl

RUN curl -L https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz | gunzip > freebayes && chmod +x freebayes

ENTRYPOINT ["/freebayes"]
