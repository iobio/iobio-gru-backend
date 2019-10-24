These instructions are for development of new backend services, using our
AWS gru backend instance.

# Preliminary

First you'll need a copy of the `dev.backend.iobio.io.pem` identity file. You
can get this from Anders.

# Launching an SSH session

```
ssh -i dev.backend.iobio.io.pem ubuntu@dev.backend.iobio.io
```

In the `ubuntu` home directory there is a `data` directory. This contains all
the necessary data (genome references, vep cache, etc) for running the backend.
This is important because gru expects `./data` to be available at runtime.
The easiest way to accomplish this is to be sure and run from the $HOME
directory.

# Create a working directory

Since this is a shared machine, you should create your own sandbax directory:

```
mkdir anders
cd anders
```

# Installing gru

First clone the repo:

```
git clone https://github.com/iobio/iobio-gru-backend
cd iobio-gru-backend
```

Then do a local install. This downloads node and all the tools:

```
make local_install
```

Finally, cd back to $HOME and run the backend (remember, you need to run
from $HOME so gru can see `$HOME/data`:

```
cd
anders/iobio-gru-backend/run_local.sh
```

You should be able to curl from your local machine now:

```
curl dev.backend.iobio.io:9001/alignmentHeader?url=http://s3.amazonaws.com/iobio/NA12878/NA12878.autsome.bam
```

You can also specify a port (default is 9001, 9001-9010 are open to the internet):

```
anders/iobio-gru-backend/run_local.sh 9002
```

# Optional - SSHFS mount

[SSHFS](https://github.com/libfuse/sshfs) is a way of mounting a remote
directory locally, using [FUSE](https://github.com/libfuse/libfuse). You could
use a normal SSH session and modify the files directly on the server using
something like VIM, but SSHFS is probably the easist method, since you can
use whatever text editor and other tools you want.

Note that anything you can do on the FUSE mount (ie git commands, copying
files, running CLI tools), you can also do with a normal SSH session. Normal
SSH will typically have much better performance. SSHFS is just for convenience
for things like text editing.

First, you'll need to install FUSE and SSHFS. You can search for how to do that
on your OS.

Once SSHFS is installed, create a directory somewhere on your local machine
where you want to mount the remote filesystem:

```
mkdir gru
```

Then, run SSHFS:

```
sshfs ubuntu@dev.backend.iobio.io: $PWD/gru -o IdentityFile=$PWD/dev.backend.iobio.io.pem
```

It's important to use an absolute path for the identity file, thus the use of
$PWD.

To umount, run `fusermount -u gru`.
