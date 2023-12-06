# Introduction

gru depends on a large data directory in order to support all functionality.
This directory contains ~128GB of data. Most of this is the VEP cache,
gnomad files, complete reference sequences, and MD5 reference sequences.

All required data is assumed to be under a `data` directory, relative to the
current working directory where gru is started from. This can be overriden with
the `--data-dir=<dir>` argument. The docker image defaults to using
`/gru_data`.

**Note:** There is some legacy cruft in the way the data volume is structured,
which can make it confusing. For example, there is `/data/data`, and also
`/data/data/references` and `/data/references`. These are all related but used
in different ways. We plan to clean this up eventually but for now you'll need
to be careful to make sure everything ends up in it's proper place.

You can use rsync to get a copy of the data directory:

## rsync

```
rsync -av rsync://data.iobio.io:9009/gru/data/gru_data_1.11.0 .
```

### Incremental updates

Note that since version `1.10.0` we're tracking changes between versions, which
allows for updates without having to download the entire directory for each
version.

For this you need 3 things:

1. The `do_update.sh` script.
2. The version you're updating from
3. The update directory that contains the changes from your current version to
   the target version.

If you're not sure what version you currently have, relatively recent versions
include a `VERSION` file you can use to determine this:

```
> cat gru_data/VERSION
> 1.10.0
```

If there's been more than one new version since the last time you updated,
you'll need to perform each update in turn. This command will list the
available updates:

```
rsync rsync://data.iobio.io:9009/gru/updates/
```

At some point it's easier to just download a fresh copy of the latest version.

Assuming you wanted to update from `1.10.0` to `1.11.0`, and you currently have
`data/gru_data_1.10.0/` in your working directory, you would do the following:

```
rsync -av rsync://data.iobio.io:9009/gru/updates/updates_1.10.0_to_1.11.0 .
rsync -av rsync://data.iobio.io:9009/gru/do_update.sh .
./do_update.sh 1.10.0 1.11.0
```

If you run into problems you can inspect `do_update.sh` to see how it works.
It's a pretty small script.

Since this uses hard links, once you've completed the update it should be safe
to delete the old version and the update directory.


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
