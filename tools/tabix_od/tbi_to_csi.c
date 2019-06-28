#include <stdio.h>
#include <stdlib.h>
#include "htslib/hts.h"

int main(int argc, char **argv)
{
    hts_idx_t *idx;
    int ret_val;

    if(argc != 3)
    {
        printf("usage: tbi_to_csi <tbi filename> <csi filename>\n");
        exit(1);
    }

    idx = hts_idx_load(argv[1], HTS_FMT_TBI);
    ret_val = hts_idx_save(idx, argv[2], HTS_FMT_CSI);

    if(ret_val < 0)
    {
        printf("Error writing %s.\n", argv[2]);
        exit(ret_val);
    }
    exit(0);
}
