#FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_109.3

USER root
WORKDIR /

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y \
        curl \
        xz-utils \
        python-pip \
        build-essential


ARG node_version=16.x
RUN curl -sL https://deb.nodesource.com/setup_${node_version} | bash - && \
    apt-get install -y nodejs

# clinReport dependencies
RUN apt-get -y install libfontconfig1

# We expect node binary to be named "node", but Ubuntu names it "nodejs", so we create a symbolic link
RUN ln -s /usr/bin/nodejs /usr/local/bin/node

#ARG gru_version=improve-logging
#RUN git clone --depth 1 --branch ${gru_version} https://github.com/iobio/iobio-gru-backend


RUN ["mkdir", "tool_bin"]

RUN curl -L https://github.com/brentp/vcfanno/releases/download/v0.3.5/vcfanno_linux64 -o tool_bin/vcfanno \
    && chmod +x tool_bin/vcfanno

RUN curl -L https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz | gunzip > tool_bin/freebayes \
        && chmod +x tool_bin/freebayes

# TODO: see if we can just copy what we need from the vep image rather than
# basing off it
#COPY --from=ensemblorg/ensembl-vep:latest /opt/vep /vep
COPY --from=quay.io/iobio/htslib:1.17 /samtools-1.17/samtools tool_bin/samtools
COPY --from=quay.io/iobio/htslib:1.17 /samtools-1.17/htslib-1.17/tabix tool_bin/tabix
COPY --from=quay.io/iobio/htslib:1.17 /samtools-1.17/htslib-1.17/bgzip tool_bin/bgzip
COPY --from=quay.io/iobio/htslib:1.17 /bcftools-1.17/bcftools tool_bin/bcftools
COPY --from=quay.io/iobio/bai-read-depther:latest /bamReadDepther tool_bin/baiReadDepther
COPY --from=quay.io/iobio/crai-read-depther:latest /craiReadDepther tool_bin/craiReadDepther
COPY --from=quay.io/iobio/bam-stats-alive:latest /bamstatsAlive/bamstatsAlive tool_bin/bamstatsAlive

COPY --from=quay.io/iobio/iobio-clin-report:latest /clin-report tool_bin/lib/clin-report
RUN cp tool_bin/lib/clin-report/clinReport.sh tool_bin/clinReport

COPY --from=quay.io/iobio/iobio-clinphen:latest /iobio-gru-backend/tool_bin/lib/clinphen/miniconda3 tool_bin/lib/clinphen/miniconda3
COPY --from=quay.io/iobio/iobio-clinphen:latest /iobio-gru-backend/tool_bin/lib/clinphen/nltk_data /usr/share/nltk_data
COPY --from=quay.io/iobio/iobio-clinphen:latest /clinphen.sh tool_bin/clinphen

COPY --from=quay.io/iobio/iobio-phenotype-extractor:latest /iobio-gru-backend/tool_bin/lib/phenotypeExtractor/miniconda3 tool_bin/lib/phenotypeExtractor/miniconda3
# TODO: Already copied from clinphen above. Maybe not safe to share the same one?
#COPY --from=quay.io/iobio/iobio-phenotype-extractor:latest /iobio-gru-backend/tool_bin/lib/phenotype-extractor/nltk_data /usr/share/nltk_data
COPY --from=quay.io/iobio/iobio-phenotype-extractor:latest /iobio-gru-backend/tool_bin/lib/phenotypeExtractor/Phenotype-extractor tool_bin/lib/phenotypeExtractor/Phenotype-extractor
COPY --from=quay.io/iobio/iobio-phenotype-extractor:latest /phenotypeExtractor.sh tool_bin/phenotypeExtractor

COPY --from=quay.io/iobio/iobio-coverage:latest /summarize_coverage.py tool_bin/summarize_coverage
COPY --from=quay.io/iobio/iobio-gene-coverage:latest /genecoverage tool_bin/geneCoverage

COPY --from=quay.io/iobio/iobio-gt-enricher:latest /usr/local/bin/gtenricher tool_bin/lib/gtEnricher/
COPY --from=quay.io/iobio/iobio-gt-enricher:htslib-1.11 /usr/local/lib/libhts.so.1.11 tool_bin/lib/gtEnricher/libhts.so.3
COPY --from=quay.io/iobio/iobio-gt-enricher:htslib-1.11 /gtEnricher.sh tool_bin/gtEnricher

COPY --from=quay.io/iobio/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/miniconda3 tool_bin/lib/knownVariants/miniconda3
COPY --from=quay.io/iobio/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/knownVariants.py tool_bin/lib/knownVariants/knownVariants.py
COPY --from=quay.io/iobio/iobio-known-variants:latest /knownVariants tool_bin/knownVariants

COPY --from=quay.io/iobio/iobio-vcf-read-depther:latest /build/vcfReadDepther tool_bin/
COPY --from=quay.io/iobio/iobio-vcf-read-depther:latest /build/lib/vcfReadDepther/vcfReadDepther tool_bin/lib/vcfReadDepther/vcfReadDepther

COPY --from=quay.io/iobio/iobio-vcf-stats-alive:latest /vcfStatsAlive tool_bin/

COPY --from=quay.io/iobio/iobio-vt:c6bd159 /vt/vt tool_bin/

EXPOSE 9001

ENTRYPOINT ["node", "/iobio-gru-backend/src/index.js"]
CMD ["--tools-dir=/iobio-gru-backend/tool_bin", "--data-dir=/gru_data"]

# These directories are required for mounting writable singularity sandbox
# images in the Utah CHPC environment, because of the directories that CHPC
# automatically binds. Hopefully this won't be necessary in the future.
RUN ["mkdir", "/uufs"]
RUN ["mkdir", "/scratch"]
RUN ["mkdir", "/ssd"]

# TODO: This is a bit hacky. Manually copying this script in order to sequence
# the docker layers such that we can rebuild code changes without
# re-downloading the static files.
COPY ./dev_tools/populate_static.sh populate_static.sh
RUN ["./populate_static.sh"]
COPY ./ /iobio-gru-backend
RUN ["mv", "/static", "/iobio-gru-backend"]
RUN ["mv", "/tool_bin", "/iobio-gru-backend"]
WORKDIR /iobio-gru-backend
RUN ["npm", "install"]
