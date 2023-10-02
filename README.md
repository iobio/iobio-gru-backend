# Run instructions

## Download a copy of the data directory

The iobio backend requires a ~128GB data directory which includes various
databases and references files.

`rsync -av rsync://data.iobio.io:9009/gru/gru_data_1.10.0/ gru_data_1.10.0`


## Run the docker image

The docker image is fairly large (~4.5GB). This creates a temporary container
with the data directory mounted

`docker run --rm -it -v $PWD/gru_data_1.10.0/:/gru_data -p 9001:9001 iobio/iobio-gru-backend:1.15.0`


## Test

In a separate terminal, the following command should return a 200 status code

`curl localhost:9001`
