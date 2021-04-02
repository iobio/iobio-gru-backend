#!/bin/bash

source /iobio-gru-backend/tool_bin/lib/phenotypeExtractor/miniconda3/bin/activate

/iobio-gru-backend/tool_bin/lib/phenotypeExtractor/Phenotype-extractor/cli.js $@
