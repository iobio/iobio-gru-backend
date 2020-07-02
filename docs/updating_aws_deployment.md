This document is specific instructions for updating the version of gru running
on the official iobio AWS deployment.

See [here](./deployment.md) for general deployment information.


# Overview 

From a high level, we have multiple gru worker nodes which handle requests.
There is an Application Load Balancer (ALB) which routes requests to the
workers. There is an autoscaling group which makes sure the proper number of
healthy nodes stays online, killing and replacing them as necessary.

In order for the autoscaling to work, we have to provide AWS with an AMI and a
launch template. When updates are made to gru, we have to update the AMI and
launch template, and then kill the old nodes (usually half at a time) so the
autoscaling replaces them with the updated version.

So the steps are:

1. Update the AMI with the latest gru

2. Update the launch template to use the new AMI

3. Kill the old nodes a few at a time so the autoscaling replaces them with the
   new ones.


# Updating the AMI

The easiest way to get a running copy of the current version of the AMI is to
simply SSH into one of the active nodes. This also ensures that you don't
pick the wrong image. Go to the EC2 instances page and pick one from there.
They are named similar to `gru-backend-worker-0.14.0`.

You'll need a copy of the `gru-backend-servers.pem` key to ssh. The command
will look something like this:

```bash
ssh -i gru-backend-servers.pem ubuntu@18.212.55.29
```

Then update gru:

```bash
cd iobio-gru-backend
git pull origin master
```

Also run the following unless you're sure you don't need to:

```bash
make clean
make
```

Scroll through the output and make sure everything downloaded correctly.

Now save a new version of the AMI from the AWS console by doing
`Actions>Image>Create Image`. Name it using the following pattern, but
increment the gru version number (and the golang and/or singularity versions if
you also update those):

```
gru-backend-0.14.0_ubuntu-18.04_golang-1.13_singularity-3.5.3
```


# Updating the launch template

Once the AMI is finished being created, use the AWS console to navigate to
the launch templates page. The launch template is named `gru-node`. Open it
and go to `Actions>Modify Template (Create new version)`.

Perform the following steps:

1. Select the `Source template version` dropdown and choose the latest version.

2. Scroll down to the AMI dropdown and select the new AMI you created above.
   Save the template version.

3. Scroll down to the resource tags section and change the `Name` tag to match
   the new version. This shouldn't affect it working but does make it easier to
   keep track of which versions we have running.

The ALB is configured to use the latest version of the launch template, so
any new nodes created after this point should use the correct AMI.


# Replacing the outdated nodes

We need to terminate the old nodes so they can be replaced with the new ones.
However, we don't want to shut them all down at once, because it takes a few
minutes for the autoscaling to detect them missing and spin up replacements.
I usually kill half of them, wait for the new ones to come up, them terminate
the other half and verify everything returns to normal.
