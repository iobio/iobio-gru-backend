#!/usr/bin/env python3

import os, argparse, threading, queue, json, math
from urllib import request
from datetime import datetime

job_queue = queue.Queue()

def worker():

    while True:
        url, parent_dir, gem_data = job_queue.get()

        print(url)

        if url.endswith('/'):
            handle_dir(url, parent_dir)
        else:
            handle_file(url, parent_dir, gem_data)

        job_queue.task_done()

def handle_dir(url, parent_dir):

    name = os.path.basename(url[:-1])
    path = os.path.join(parent_dir, name)

    try:
        os.makedirs(path)
    except:
        pass

    gem_url = url + '.gemdrive-ls.json'

    res = request.urlopen(gem_url)
    body = res.read()
    gem_dir = json.loads(body)

    if 'children' not in gem_dir:
        return

    for child_name in gem_dir['children']:
        child = gem_dir['children'][child_name]
        child_url = url + child_name
        job_queue.put((child_url, path, child))

def handle_file(url, parent_dir, gem_data):

    name = os.path.basename(url)
    path = os.path.join(parent_dir, name)

    try:
        stat = os.stat(path)
    except:
        stat = None

    size = 0
    mod_time = ''
    if stat:
        size = stat.st_size
        mod_time = datetime.utcfromtimestamp(stat.st_mtime).replace(microsecond=0).isoformat() + 'Z'

    utc_dt = datetime.strptime(gem_data['modTime'], '%Y-%m-%dT%H:%M:%SZ')
    mtime = math.floor((utc_dt - datetime(1970, 1, 1)).total_seconds())

    needs_update = size != gem_data['size'] or mod_time != gem_data['modTime']

    if needs_update:
        res = request.urlopen(url)
        with open(path, 'wb') as f:
            while chunk := res.read(4096):
                f.write(chunk)
        stat = os.stat(path)
        os.utime(path, (stat.st_atime, mtime))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('url', help='GemDrive directory URL')
    parser.add_argument('--num-workers', type=int, help='Number of worker threads', default=4)
    parser.add_argument('--out-dir', help='Output directory', default=os.getcwd())
    args = parser.parse_args()

    for w in range(args.num_workers):
        threading.Thread(target=worker, daemon=True).start()

    job_queue.put((args.url, args.out_dir, {}))
    job_queue.join()
