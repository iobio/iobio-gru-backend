# Introduction

gru depends on a large data directory in order to support all functionality.
This directory contains ~120GB of data. Most of this is the VEP cache,
complete reference sequences, and MD5 reference sequences.

All required data is assumed to be under a `data` directory, relative to the
current working directory where gru is started from (usually with
`run_local.sh`). This can be overriden with the `--out-dir=<dir>` argument.
(note that we only recently fixed some bugs with `--out-dir`, and it's not
fully tested yet. Please let us know if you run into issues with this
functionality).

**Note:** There is some legacy cruft in the way the data volume is structured,
which can make it confusing. For example, there is `/data/data`, and also
`/data/data/references` and `/data/references`. These are all related but used
in different ways. We plan to clean this up eventually but for now you'll need
to be careful to make sure everything ends up in it's proper place.

There are 2 simple ways to get a copy of the data directory:


## AWS snapshot

If you're on AWS, you can create a volume from our public snapshot. The current
gru version is 0.28.0, which is snapshot ID `snap-02eb780cd49ee41fe`.


## GemDrive

Download a copy from our GemDrive instance. [GemDrive] is a simple HTTP
protocol that allows for recursively downloading directories.

You can explore the data from a web browser here:

https://gemdrive.io/apps/delver/?drive=https://gemdrive.iobio.io&path=/

You can download the files from that interface (for example navigate to
`/gru-0.28.0/iobio-gru-backend/data` to see the data directory for the 0.28.0
release), but it's not efficient for datasets this large. We have a python3
script in `dev_tools/download.py` that can recursively download the data
directory using the following command:

```bash
python3 dev_tools/download.py
```

There are a few arguments to override defaults:


```bash
python3 dev_tools/download.py --gru-version 0.27.0 --out-dir out --dry-run
```


# Manual downloads

Two of the biggest pieces of the data directory are the VEP cache and the
MD5 reference cache, both of which are available publicly. If you can download
them manually from their source locations, it reduces the load on our servers.
Instructions are below.

Note that you can populate these directories manually, then run the download.py
script with the same output `data` directory. As long as the file modification
times match, it will only overwrite missing data from the GemDrive.

## VEP cache

gru requires a vep cache installed under `data/vep-cache`. You can follow the
instructions [here][0] to manually download the appropriate cache data. The
version we're currently using is **101**. You will need the following files
(note that these are the indexed versions):

* ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_refseq_vep_101_GRCh37.tar.gz
* ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_refseq_vep_101_GRCh38.tar.gz
* ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_vep_101_GRCh37.tar.gz
* ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_vep_101_GRCh38.tar.gz

These need to be decompressed into the following directory structure:

```
data/
  vep-cache/
    homo_sapiens/
      101_GRCh37/
      101_GRCh38/
    homo_sapiens_refseq/
      101_GRCh37/
      101_GRCh38/
```

## MD5 Reference Cache

gru uses reference sequences (in several different formats) for many endpoints.
The most obvious is that an MD5 reference cache is required to operate on CRAM
files. By default, if an MD5 reference isn't available locally, samtools will
reach out to the EBI servers to attempt to download the reference at runtime.
This can be really slow. Fortunately, the cache can be pre-downloaded, using
the `seq_cache_populate.pl` [perl script][1]. The script runs on a fasta
reference file. It goes through the file and creates an MD5 of each entry,
then copies the contents of the entry into a special directory under that
filename. See [here][2] for more details.

If your reference sequence is named `human_g1k_v38_decoy_phix.fasta`, you
would generate the MD5 cache for GRCh38 like this:

`./seq_cache_populate.pl human_g1k_v38_decoy_phix.fasta -root md5_reference_cache`

Since everything is based off hashes, you can use the same output directory
(-root) for builds 37 and 38). Move that directory under your data directory
like this:

```
data/
  md5_reference_cache/
```

[0]: https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

[1]: https://github.com/samtools/samtools/blob/develop/misc/seq_cache_populate.pl

[2]: ./handling_cram_references.md

[GemDrive]: https://github.com/gemdrive
