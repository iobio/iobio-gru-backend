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


## Data directory

gru expects a copy of the data directory to exist at `./data`, relative to
the working directory where you run gru.

See [here](populating_data_directory.md) for more information on creating the
data directory.


## Running

You can run `run_local.sh` to start the server. Make sure the data directory
exists in the working directory you run from. The default port is 9001. You
can override it like this:

```
./run_local.sh --port=9002`
```


# Production Deployment

For a production deployment, you'll want to have a few extra components:

1. Multiple instances

2. A load balancer to hand out requests to instances

3. Monitoring to alert/shutdown/replace failed instances

4. (optional) auto scaling to add/remove instances based on system load


## iobio AWS deployment

The public production deployment of gru runs on Amazon AWS. This section
details some of the setup we use.


### Data directory EBS

We maintain an EBS volume which contains a copy of the latest gru, including
the code (this repository) located in `iobio-gru-backend/`, and the data
directory in `iobio-gru-backend/data/`.

The current version of the EBS is 0.28.0. You can obtain a copy by using the
following AWS snapshot:

`snap-02eb780cd49ee41fe`

Let's say you create a copy of that EBS and mount it at `/mnt/gru`. You can
then run gru like this:

```
cd /mnt/gru/iobio-gru-backend
./run_local.sh --port=9002
```

Note that EBS snapshots are stored on Amazon S3, which is much slower than EBS.
When you create a new EBS volume from a snapshot, files are lazy-loaded from S3
on demand. This can cause serious performance issues the first time each file
is accessed. This means you probably don't want to recreate gru EBS volumes
from scratch too often if you can avoid it. See [here][0].


### VM Images

We maintain an AMI for launching gru instances. The current AMI ID is:

```
ami-08333660d9a8a8598
```

If you want to use that AMI, note that it's decoupled from any specific version
of the EBS volume, so you'll need to make an AWS [launch template][1] to launch
the image with a copy of the snapshot version you want attached. I didn't see a
way to share our launch template publicly, but it's not doing much so hopefully
won't be too difficult to make your own.


#### How we make the AMI

First, we follow all the steps above for a fresh Ubuntu 18.04 VM, including
enabling the systemd service.

Then, an fstab entry is created for the data EBS, like so:

```
# /etc/fstab
UUID=946d678a-9da6-41b1-9564-26c8fde76852 /home/ubuntu/data ext4 ro,suid,dev,auto 0 0
```

You can find the UUID by running `sudo blkid` after attaching EBS volume.

From that state, we create an AMI which can be used to start fresh gru
instances.


### Load balancer

We set up an application load balancer, with health checks going to the root
HTTP path of the instances.


### Auto scaling group

We use an auto scaling group, in tandem with the load balancer. This requires
creating a launch template mentioned above.

Since the systemd starts gru on boot, the auto scaling is able to replace
unhealthy nodes easily. It's also easy to manually delete a node and it will
be replaced automatically.


[0]: https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-initialize.html

[1]: https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-launch-templates.html
