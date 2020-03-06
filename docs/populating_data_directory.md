gru depends on a large data directory in order to support all functionality.
This directory contains ~100GB of data. Most of this is the VEP cache,
complete reference sequences, and MD5 reference sequences.

All required data is assumed to be under a `data` directory, relative to the
current working directory.


# VEP cache

gru requires a vep cache installed under `data/vep-cache`. You can follow the
instructions [here][0] to manually download the appropriate cache data. The
version we're currently using is **97**. You will need the following files
(note that these are the indexed versions):

* ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_refseq_vep_97_GRCh37.tar.gz
* ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_refseq_vep_97_GRCh38.tar.gz
* ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_97_GRCh37.tar.gz
* ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_97_GRCh38.tar.gz

These need to be decompressed into the following directory structure:

```
data/
  vep-cache/
    homo_sapiens/
      97_GRCh37/
      97_GRCh38/
    homo_sapiens_refseq/
      97_GRCh37/
      97_GRCh38/
```

[0]: https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
