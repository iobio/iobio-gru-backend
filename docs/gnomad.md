This will tell you how to create the gnomad VCF files used by iobio


# Download gnomad

Get the full gnomad files for each chromosome from
<https://gnomad.broadinstitute.org/downloads>.


# Filter files

The gnomad files are large, but we don't need all the annotations. So we need
to filter out the unused data.

You first need to do the equivalent of `dev_tools/filter_gnomad_vcfs.sh`
followed by `dev_tools/merge_gnomad_vcfs.sh` to create the iobio gnomad
VCF and index files.

`dev_tools/filter_gnomad_vcfs.sh` is written so it can be run on the CHPC slurm
cluster. It should be straight forward to adapt it to work on a different slurm
cluster, or in a standard terminal session, though the latter would probably
take a long time.
