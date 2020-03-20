gru depends on a large data directory in order to support all functionality.
This directory contains ~100GB of data. Most of this is the VEP cache,
complete reference sequences, and MD5 reference sequences.

All required data is assumed to be under a `data` directory, relative to the
current working directory where gru is started from (usually with
`run_local.sh`).

**Note:** There is some legacy cruft in the way the data volume is structured,
which can make it confusing. For example, there is `/data/data`, and also
`/data/data/references` and `/data/references`. These are all related but used
in different ways. We plan to clean this up eventually but for now you'll need
to be careful to make sure everything ends up in it's proper place.


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

**Note:** The files don't unzip to the proper directory names. You'll need to
rename them to match the names above.


# MD5 Reference Cache

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

If you don't have a reference, one is included in the Google Drive data
folder.

# Other files

gru requires several other directories and files to operate. These can be
downloaded from the Google Drive share. The final top-level data directory
should look like this:


```
data/
  data/
  gene2pheno/
  geneinfo/
  genomebuild/
  gnomad_header.txt
  md5_reference_cache/
  references/
  vep-cache/
```

[0]: https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

[1]: https://github.com/samtools/samtools/blob/develop/misc/seq_cache_populate.pl

[2]: ./handling_cram_references.md
