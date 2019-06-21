#!/bin/bash

workers=$(./get_aws_addresses.py)

for worker in ${workers}
do
    echo $worker

    echo Copying files to ${worker}
    rsync -av -e "ssh -i ${SSH_KEY_FILE} -oStrictHostKeyChecking=no" --exclude=".git" --delete ./* ubuntu@${worker}:iobio-backend
    
    echo Starting server
    ssh -i ${SSH_KEY_FILE} ubuntu@${worker} killall node
    ssh -i ${SSH_KEY_FILE} ubuntu@${worker} 'cd iobio-backend; PATH=./tools:$PATH nohup node/bin/node index.js > log.out 2> log.err < /dev/null &'
    
    echo Done

done
