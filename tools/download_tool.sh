#!/bin/bash
tool=${1}

wget https://iobio.s3.amazonaws.com/backend_tools/${tool} -O tool_bin/${tool}
chmod +x tool_bin/${tool}
