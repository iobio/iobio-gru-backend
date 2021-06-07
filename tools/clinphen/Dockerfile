FROM ubuntu:18.04

RUN apt-get update && apt-get install -y curl

RUN mkdir -p /iobio-gru-backend/tool_bin/lib/clinphen

WORKDIR /iobio-gru-backend/tool_bin/lib/clinphen

COPY build.sh /
RUN chmod +x /build.sh
RUN /build.sh

COPY clinphen.sh /
RUN chmod +x /clinphen.sh
