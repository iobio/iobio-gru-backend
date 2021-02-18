#!/bin/bash

DIR=/iobio-gru-backend/tool_bin/lib/clin-report

export PATH=$DIR/node/bin:$PATH

node $DIR/index.js $@
