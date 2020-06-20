# Introduction

gru requires 2 primary pieces to run:

1. A node.js server

2. A data directory (~100GB)

Simple deployments will only require a single instance, but the system is
designed to be scaled horizontally as needed. Every endpoint is stateless. The
trickiest part of deployment is that every instance needs access to a copy of
the data directory. This can be accomplished by creating a separate copy of the
data directory for each instance, or using shared filesystems like NFS, Lustre,
CEPH, etc.


# Prerequisites

The host machine must run linux with singularity installed. The simplest way
to get set up is to clone the gru repo, then use the `install_go.sh` script
to install golang, which is required to build singularity. Then run
`install_singularity.sh` to complete setup.

## Steps to set up fresh Ubuntu 18.04 instance:

1. `git clone https://github.com/iobio/iobio-gru-backend`

2. `iobio-gru-backend/install_go.sh`

3. `sudo apt-get install build-essential libssl-dev uuid-dev`

4. `iobio-gru-backend/install_singularity.sh`

5. Optionally, you can create a systemd service to run gru automatically on
   boot, and to restart it if it crashes. A systemd service file is included in
   the repo at `systemd/iobio-gru-backend.service`. Copy it to
   `/etc/systemd/system` and run `sudo systemctl enable iobio-gru-backend`
   (this sets it up to start on boot), then run `sudo systemctl start
   iobio-gru-backend`.


# Simple deployment


## Server

In the simplest case, all that's necessary for setting up gru is to clone the
repo, then run `make`. If something goes wrong, it should be fairly easy to
manually run the steps in the Makefile. There's not much going on in there.

Note that a lot of bioinformatics CLI tools are required by gru at runtime.
These tools must be in `tool_bin`.  By default, we pre-build these tools and
store them on a public S3 bucket, and the Makefile downloads them.

Recipes to build the tools from scratch are included in the `tools/`
directory. Most of them are packaged as singularity containers. As an example,
to build samtools you would run the following from `tools/samtools/`:

```
sudo singularity build samtools samtools.def
```

And then copy the samtools image to `tool_bin/`. It can mostly be treated as
a static binary, but whatever system you run it on must have singularity
installed.


## Data directory

gru expects a copy of the data directory to exist at `./data`, relative to
the working directory.

See [here](populating_data_directory.md) for more information on creating the
data directory.


## Running

You can run `run_local.sh` to start the server. Make sure the data directory
exists in the working directory you run from.

You can also run `src/index.js` manually. The key thing is properly setting up
your `PATH`. See `run_local.sh` for details.


# Production Deployment

For a production deployment, you'll want to have a few extra components:

1. Multiple instances

2. A load balancer to hand out requests to instances

3. Monitoring to alert/shutdown/replace failed instances

4. (optional) auto scaling to add/remove instances based on system load


## iobio AWS deployment

The main production deployment of gru runs on Amazon AWS. This section
details some of the setup we use.


### Data directory EBS

We maintain an EBS volume which contains a copy of the latest data
directory. Whenever we update the data directory, we create a snapshot of
that EBS.


### VM Images

We maintain an AMI for launching gru instances.

First, we follow all the steps above for a fresh Ubuntu 18.04 VM, including
enabling the systemd service.

Then, an fstab entry is created for the data EBS, like so:

```
# /etc/fstab
UUID=946d678a-9da6-41b1-9564-26c8fde76852 /home/ubuntu/data ext4 ro,suid,dev,auto 0 0
```

You can find the UUID by running `sudo blkid` after attaching the EBS volume.


From that state, we create an AMI which can be used to start fresh gru
instances.

Updated AMIs can be created by modifying existing VMs and creating a new AMI,
or by running these steps from scratch. The `deploy_aws.sh` script is useful
if you need to update a VM (or multiple VMs) while it's running, as long as
you don't need to update the data volume at the same time. You'll want to
read the script to see how it works, here's how we run it as an example:

```
SSH_KEY_FILE=$HOME/path/to/sshkey.pem ./deploy_aws.sh prod
```

### Load balancer

We set up an application load balancer, with health checks going to the root
HTTP path of the instances.


### Auto scaling group

We use an auto scaling group, in tandem with the load balancer. This requires
creating a launch template to specify things like the AMI and instance types
to use.

Since the systemd starts gru on boot, the auto scaling is able to replace
unhealthy nodes easily. It's also easy to manually delete a node and it will
be replaced automatically.
