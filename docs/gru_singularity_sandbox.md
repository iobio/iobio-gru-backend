# Creating the sandbox

You'll need a `gru.sif` singularity image to get started. Once you have that,
you can run the following commands.

First load the singularity module (if running in CHPC):

```
module load singularity
```

Then create the sandbox directory:


```
singularity build --sandbox gru_sandbox/ gru.sif
```



# Running gru

The command to run gru is a bit complicated, but the key ingredients are:

1. Bind a data directory to the container.
2. Run gru on the desired port
3. Tell gru where to find the data directory

Here's what it ends up looking like:

```
singularity exec --bind /ssd/gru/gru-1.0.0/data/:/data gru.sif node /iobio-gru-backend/src/index.js --port=9003 --data-dir=/data
```


# Changing the code

The easiest way to make code modifications is to change the files directly
in the sandbox directory, rather than shelling in to the sandbox. This lets
you use your editor of choice and your current git config.

Then you can go on GitHub and create a pull request from your branch.
Currently you'll want to submit the PR to
`https://github.com/anderspitman/iobio-gru-backend`.
