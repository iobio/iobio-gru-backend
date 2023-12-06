# Run instructions

## Download a copy of the data directory

The iobio backend requires a ~128GB data directory which includes various
databases and references files.

```
rsync -av rsync://data.iobio.io:9009/gru/data/gru_data_1.11.0 .
```

Note that since `gru_data_1.10.0` it's possible to upgrade data directory
versions without a complete download. See
[here](docs/populating_data_directory.md#incremental-updates) for details.


## Run the docker image

The docker image is fairly large (~4.5GB). This creates a temporary container
with the data directory mounted

`docker run --rm -it -v $PWD/gru_data_1.11.0/:/gru_data -p 9001:9001 iobio/iobio-gru-backend:1.17.0`


## Test

In a separate terminal, the following command should return a 200 status code

`curl localhost:9001`
