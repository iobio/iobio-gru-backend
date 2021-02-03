FROM centos:7

#RUN apt-get update && apt-get -y install curl autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

RUN yum update -y && yum install -y autoconf automake make gcc zlib-devel bzip2 bzip2-devel xz-devel curl-devel

ARG htslib_version=1.11

COPY build.sh /

RUN chmod +x build.sh
RUN /build.sh ${htslib_version}
