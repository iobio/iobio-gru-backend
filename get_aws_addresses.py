#!/usr/bin/env python
import boto3
from pprint import pprint
import sys

target = sys.argv[1]

target_name = 'gru-backend-worker-0.10.0'
if target == 'prod':
    target_name = 'gru-backend-worker-0.10.0'

client = boto3.client('ec2')
response = client.describe_instances()
for reservation in response['Reservations']:
    for instance in reservation['Instances']:
        if 'Tags' in instance:
            for tag in instance['Tags']:
                if instance['State']['Name'] == 'running' and tag['Key'] == 'Name' and tag['Value'] == target_name:
                    print(instance['PublicDnsName'])
