#FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_101.0

USER root
WORKDIR /

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y \
        curl \
        git \
        xz-utils \
        python-pip \
        build-essential


ARG node_version=12.x
RUN curl -sL https://deb.nodesource.com/setup_${node_version} | bash - && \
    apt-get install -y nodejs

# clinReport dependencies
RUN apt-get -y install libfontconfig1

# We expect node binary to be named "node", but Ubuntu names it "nodejs", so we create a symbolic link
RUN ln -s /usr/bin/nodejs /usr/local/bin/node

ARG gru_version=dev
RUN git clone --depth 1 --branch ${gru_version} https://github.com/anderspitman/iobio-gru-backend

WORKDIR /iobio-gru-backend

RUN ["./dev_tools/populate_static.sh"]

RUN ["npm", "install"]

RUN ["mkdir", "tool_bin"]

# TODO: see if we can just copy what we need from the vep image rather than
# basing off it
#COPY --from=ensemblorg/ensembl-vep:latest /opt/vep /vep
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/samtools tool_bin/samtools-1.11
RUN ln -s $PWD/tool_bin/samtools-1.11 $PWD/tool_bin/samtools
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/tabix tool_bin/tabix-1.11
RUN ln -s $PWD/tool_bin/tabix-1.11 $PWD/tool_bin/tabix
COPY --from=anderspitman/htslib:1.11 /samtools-1.11/htslib-1.11/bgzip tool_bin/bgzip-1.11
RUN ln -s $PWD/tool_bin/bgzip-1.11 $PWD/tool_bin/bgzip
COPY --from=anderspitman/htslib:1.11 /bcftools-1.11/bcftools tool_bin/bcftools-1.11
RUN ln -s $PWD/tool_bin/bcftools-1.11 $PWD/tool_bin/bcftools
COPY --from=anderspitman/bai-read-depther:latest /bamReadDepther tool_bin/baiReadDepther
COPY --from=anderspitman/crai-read-depther:latest /craiReadDepther tool_bin/craiReadDepther
COPY --from=anderspitman/bam-stats-alive:latest /bamstatsAlive/bamstatsAlive tool_bin/bamstatsAlive

COPY --from=anderspitman/iobio-clin-report:latest /clin-report tool_bin/lib/clin-report
RUN cp tool_bin/lib/clin-report/clinReport.sh tool_bin/clinReport

COPY --from=anderspitman/iobio-clinphen:latest /iobio-gru-backend/tool_bin/lib/clinphen/miniconda3 tool_bin/lib/clinphen/miniconda3
COPY --from=anderspitman/iobio-clinphen:latest /iobio-gru-backend/tool_bin/lib/clinphen/nltk_data /usr/share/nltk_data
COPY --from=anderspitman/iobio-clinphen:latest /clinphen.sh tool_bin/clinphen

COPY --from=anderspitman/iobio-coverage:latest /summarize_coverage.py tool_bin/summarize_coverage
COPY --from=anderspitman/iobio-gene-coverage:latest /genecoverage tool_bin/geneCoverage

COPY --from=anderspitman/iobio-gt-enricher:latest /usr/local/bin/gtenricher tool_bin/lib/gtEnricher/
COPY --from=anderspitman/iobio-gt-enricher:latest /usr/local/lib/libhts.so.1.11 tool_bin/lib/gtEnricher/libhts.so.3
COPY --from=anderspitman/iobio-gt-enricher:latest /gtEnricher.sh tool_bin/gtEnricher

COPY --from=anderspitman/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/miniconda3 tool_bin/lib/knownVariants/miniconda3
COPY --from=anderspitman/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/knownVariants.py tool_bin/lib/knownVariants/knownVariants.py
COPY --from=anderspitman/iobio-known-variants:latest /knownVariants tool_bin/knownVariants

COPY --from=anderspitman/iobio-vcf-read-depther:latest /build/vcfReadDepther tool_bin/
COPY --from=anderspitman/iobio-vcf-read-depther:latest /build/lib/vcfReadDepther/vcfReadDepther tool_bin/lib/vcfReadDepther/vcfReadDepther

COPY --from=anderspitman/iobio-vcf-stats-alive:latest /vcfStatsAlive tool_bin/

COPY --from=anderspitman/iobio-vt:latest /vt/vt tool_bin/

COPY --from=anderspitman/iobio-freebayes:latest /freebayes tool_bin/

EXPOSE 9001

#ENTRYPOINT ["node", "/iobio-gru-backend/src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin"]
CMD ["node", "/iobio-gru-backend/src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin"]
