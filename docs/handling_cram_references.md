CRAM files require a reference sequence for reading and writing. This document
is intended as a central knowledge base about how this works with samtools.

When samtools is provided a CRAM file, by default samtools will search for
reference sequences in the following order (see
[the HTSLIB docs](http://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES)):


1. Use any local file specified by the command line options (eg -T).
2. Look for MD5 via REF_CACHE environment variable.
3. Look for MD5 in each element of the REF_PATH environment variable.
4. Look for a local file listed in the UR: header tag. 

For an explanation for REF_CACHE and REF_PATH, see
[here](http://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES) and
[here](http://www.htslib.org/workflow/#mapping_to_cram).


The default location for the cache on Linux is `$HOME/.cache/hts-ref` (see
[this mailing list message](https://sourceforge.net/p/samtools/mailman/message/34781943/)).

In the default configuration, if samtools is run on a CRAM file and no cache
exists (and no reference is provided on the command line), samtools will
attempt to automatically use the MD5s from the CRAM header to download the
proper reference sequences from EBI. For example, the reference sequence
for MD5 `ff2e5f34be64b46f1e15d1091ea54645` can be downloaded manually like this:

```bash
curl https://www.ebi.ac.uk/ena/cram/md5/ff2e5f34be64b46f1e15d1091ea54645
```

samtools will store the references in `$HOME/.cache/hts-ref`, using a format
that avoids having too many files in a single directory. The format is to use
the first 2 characters of the MD5 as the first directory, and the 3rd and 4th
characters as the second directory, and the remainder as the filename. So the
example reference above with MD5 `ff2e5f34be64b46f1e15d1091ea54645` would be
stored in `$HOME/.cache/hts-ref/ff/2e/5f34be64b46f1e15d1091ea54645`.

On subsequent runs, the cache directory is checked before attempting to
retrieve the reference over the internet.



# Reference links

* [How the MD5 is created (see section 1.3.1)](https://samtools.github.io/hts-specs/SAMv1.pdf)
