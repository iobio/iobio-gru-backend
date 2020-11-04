#!/bin/bash

trap 'kill $(jobs -p)' EXIT

# Incantation taken from https://stackoverflow.com/a/246128/943814
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

$SCRIPT_DIR/node/bin/node $SCRIPT_DIR/src/index.js $@
