#FROM ubuntu:18.04
FROM ensemblorg/ensembl-vep:release_105.0

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

COPY ./ /iobio-gru-backend

WORKDIR /iobio-gru-backend

RUN ["./dev_tools/populate_static.sh"]

RUN ["npm", "install"]

RUN ["mkdir", "tool_bin"]

RUN curl -L https://github.com/brentp/vcfanno/releases/download/v0.3.3/vcfanno_linux64 -o tool_bin/vcfanno \
    && chmod +x tool_bin/vcfanno

# TODO: see if we can just copy what we need from the vep image rather than
# basing off it
#COPY --from=ensemblorg/ensembl-vep:latest /opt/vep /vep
COPY --from=quay.io/iobio/htslib:1.11 /samtools-1.11/samtools tool_bin/samtools-1.11
RUN ln -s $PWD/tool_bin/samtools-1.11 $PWD/tool_bin/samtools
COPY --from=quay.io/iobio/htslib:1.11 /samtools-1.11/htslib-1.11/tabix tool_bin/tabix-1.11
RUN ln -s $PWD/tool_bin/tabix-1.11 $PWD/tool_bin/tabix
COPY --from=quay.io/iobio/htslib:1.11 /samtools-1.11/htslib-1.11/bgzip tool_bin/bgzip-1.11
RUN ln -s $PWD/tool_bin/bgzip-1.11 $PWD/tool_bin/bgzip
COPY --from=quay.io/iobio/htslib:1.11 /bcftools-1.11/bcftools tool_bin/bcftools-1.11
RUN ln -s $PWD/tool_bin/bcftools-1.11 $PWD/tool_bin/bcftools
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
COPY --from=quay.io/iobio/iobio-gt-enricher:latest /usr/local/lib/libhts.so.1.11 tool_bin/lib/gtEnricher/libhts.so.3
COPY --from=quay.io/iobio/iobio-gt-enricher:latest /gtEnricher.sh tool_bin/gtEnricher

COPY --from=quay.io/iobio/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/miniconda3 tool_bin/lib/knownVariants/miniconda3
COPY --from=quay.io/iobio/iobio-known-variants:latest /iobio-gru-backend/tool_bin/lib/knownVariants/knownVariants.py tool_bin/lib/knownVariants/knownVariants.py
COPY --from=quay.io/iobio/iobio-known-variants:latest /knownVariants tool_bin/knownVariants

COPY --from=quay.io/iobio/iobio-vcf-read-depther:latest /build/vcfReadDepther tool_bin/
COPY --from=quay.io/iobio/iobio-vcf-read-depther:latest /build/lib/vcfReadDepther/vcfReadDepther tool_bin/lib/vcfReadDepther/vcfReadDepther

COPY --from=quay.io/iobio/iobio-vcf-stats-alive:latest /vcfStatsAlive tool_bin/

COPY --from=quay.io/iobio/iobio-vt:latest /vt/vt tool_bin/

COPY --from=quay.io/iobio/iobio-freebayes:latest /freebayes tool_bin/

EXPOSE 9001

#ENTRYPOINT ["node", "/iobio-gru-backend/src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin"]
CMD ["node", "/iobio-gru-backend/src/index.js", "--tools-dir=/iobio-gru-backend/tool_bin", "--data-dir=/data"]

# These directories are required for mounting writable singularity sandbox
# images in the Utah CHPC environment, because of the directories that CHPC
# automatically binds. Hopefully this won't be necessary in the future.
RUN ["mkdir", "/uufs"]
RUN ["mkdir", "/scratch"]
RUN ["mkdir", "/ssd"]
