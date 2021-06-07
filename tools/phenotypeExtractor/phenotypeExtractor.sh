#!/bin/bash

#source /iobio-gru-backend/tool_bin/lib/phenotypeExtractor/miniconda3/bin/activate

export PATH=/iobio-gru-backend/tool_bin/lib/phenotypeExtractor/miniconda3/bin:$PATH

/iobio-gru-backend/tool_bin/lib/phenotypeExtractor/Phenotype-extractor/cli.js $@
