FROM ubuntu:18.04

RUN apt-get update && apt-get install -y curl

RUN mkdir -p /iobio-gru-backend/tool_bin/lib/knownVariants

WORKDIR /iobio-gru-backend/tool_bin/lib/knownVariants

COPY knownVariants.py ./

COPY build.sh /
RUN chmod +x /build.sh
RUN /build.sh

COPY knownVariants.sh /knownVariants
RUN chmod +x /knownVariants
