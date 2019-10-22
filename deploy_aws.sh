#!/bin/bash

# Example usage:
# SSH_KEY_FILE=iobioServers.cer ./deploy_aws.sh

SSH_OPTIONS="-oStrictHostKeyChecking=no"

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
    workers=$(./get_aws_addresses.py)
fi

for worker in ${workers}
do
    echo $worker

    echo Copying files to ${worker}
    rsync -av -e "ssh ${SSH_OPTIONS}" \
        --exclude=".git" --exclude="tools" --exclude="data" --exclude="vep-cache" \
        --delete ./* ubuntu@${worker}:iobio-backend
    
    echo Starting server
    ssh $SSH_OPTIONS ubuntu@$worker "killall node"
    ssh $SSH_OPTIONS ubuntu@$worker 'PATH=$PWD/iobio-backend/tool_bin:$PATH nohup $PWD/iobio-backend/node/bin/node $PWD/iobio-backend/src/index.js > log.out 2> log.err < /dev/null &'
    #ssh $SSH_OPTIONS ubuntu@$worker "ln -s -f /mnt/data/data iobio-backend/"
    #ssh $SSH_OPTIONS ubuntu@$worker "ln -s -f /mnt/data/vep-cache iobio-backend/"

    printf "Done\n\n"
done
