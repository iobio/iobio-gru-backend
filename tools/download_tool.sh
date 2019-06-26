#!/bin/bash
tool=${1}

curl https://iobio.s3.amazonaws.com/backend_tools/${tool} > tool_bin/${tool}
chmod +x tool_bin/${tool}
