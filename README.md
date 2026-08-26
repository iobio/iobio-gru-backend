# Run instructions

## Quickstart

Install [Docker](https://docs.docker.com/engine/install/) and [rclone](https://rclone.org/install/). GRU requires a roughly 128 GB data directory containing databases and reference files. Download data version `2.0.0` from `files.iobio.io`:

```bash
rclone sync --progress --http-url https://files.iobio.io \
  :http:gru_data/data/gru_data_2.0.0/ \
  gru_data_2.0.0
```

Running the command again synchronizes an existing local copy rather than downloading unchanged files again.

Start GRU with the downloaded data mounted:

```bash
docker run --rm -it \
  -v "$PWD/gru_data_2.0.0:/gru_data" \
  -p 9001:9001 \
  docker.io/iobio/iobio-gru-backend:2.0.0
```

In another terminal, verify that GRU is running:

```bash
curl http://localhost:9001
```

## Runtime configuration

GRU reads `config.json`, applies environment-variable overrides, and serves the effective configuration to hosted applications at `/config.json`. Configuration keys use snake_case. Route mounts are configured with `path_prefix`, for example:

```json
{
  "gene": {
    "path_prefix": "/gene",
    "default_mode": "advanced"
  },
  "backend": {
    "path_prefix": "/"
  }
}
```

Environment variables beginning with `IOBIO_GENE_`, `IOBIO_BAM_`, or `IOBIO_BACKEND_` are published under the corresponding configuration section. The suffix is lowercased without otherwise changing its snake_case form:

```text
IOBIO_GENE_PATH_PREFIX=/gene       -> gene.path_prefix
IOBIO_GENE_DEFAULT_MODE=simple     -> gene.default_mode
IOBIO_GENE_SHOW_INTRO=true         -> gene.show_intro
IOBIO_BACKEND_ORIGIN=https://...   -> backend.origin
```

Exact values `true` and `false` become JSON booleans. All other environment values remain strings. These namespaces are public application configuration because the resulting values are returned by `/config.json`.
