#!/usr/bin/env python
import boto3
from pprint import pprint

client = boto3.client('ec2')
response = client.describe_instances()
for reservation in response['Reservations']:
    for instance in reservation['Instances']:
        for tag in instance['Tags']:
            if instance['State']['Name'] == 'running' and tag['Key'] == 'Name' and tag['Value'].startswith('backend-worker'):
                print(instance['PublicDnsName'])
