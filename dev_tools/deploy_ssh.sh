#!/bin/bash

# Example usage:
# SSH_KEY_FILE=iobioServers.cer ./deploy_ssh.sh stage

target=$1

SSH_OPTIONS="-oStrictHostKeyChecking=no"

if [ "$target" != "stage" ] && [ "$target" != "prod" ]
then
  echo "Usage: deploy_aws.sh [stage|prod]"
  exit 1
fi

if [ -n "${SSH_KEY_FILE}" ]
then
    SSH_OPTIONS="${SSH_OPTIONS} -i ${SSH_KEY_FILE}"
    echo $SSH_OPTIONS
fi

echo "SSH options:"
echo ${SSH_OPTIONS}


if [ -n "${NODES}" ]
then
    workers=${NODES}
else
    workers=$(./get_aws_addresses.py ${target})
fi

for worker in ${workers}
do
    echo $worker

    echo Starting server
    ssh $SSH_OPTIONS ubuntu@$worker "killall run_local.sh"
    ssh $SSH_OPTIONS ubuntu@$worker '$PWD/iobio-backend/run_local.sh > log.out 2> log.err < /dev/null &'

    printf "Done\n\n"
done
