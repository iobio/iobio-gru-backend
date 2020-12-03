#!/usr/bin/env python3

import os, argparse, threading, queue, json, math
from urllib import request
from datetime import datetime

job_queue = queue.Queue(maxsize=8)

def traverse(url, parent_dir):

    print(url)

    name = os.path.basename(url[:-1])
    path = os.path.join(parent_dir, name)

    try:
        os.makedirs(path)
    except:
        pass

    gem_url = url + 'gemdrive/meta.json'

    res = request.urlopen(gem_url)
    body = res.read()
    gem_dir = json.loads(body)

    if 'children' not in gem_dir:
        return

    for child_name in gem_dir['children']:
        child = gem_dir['children'][child_name]
        child_url = url + child_name
        if child_url.endswith('/'):
            traverse(child_url, path)
        else:
            job_queue.put((child_url, path, child))


def downloader():
    while True:
        url, parent_dir, gem_data = job_queue.get()

        print(url)

        handle_file(url, parent_dir, gem_data)

        job_queue.task_done()


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
            while True:
                chunk = res.read(4096)
                if not chunk:
                    break
                f.write(chunk)
        stat = os.stat(path)

        if stat.st_size != gem_data['size']:
            print("Sizes don't match", url)

        os.utime(path, (stat.st_atime, mtime))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('url', help='GemDrive directory URL')
    parser.add_argument('--num-workers', type=int, help='Number of worker threads', default=4)
    parser.add_argument('--out-dir', help='Output directory', default=os.getcwd())
    args = parser.parse_args()

    for w in range(args.num_workers):
        threading.Thread(target=downloader, daemon=True).start()

    traverse(args.url, args.out_dir)

    job_queue.join()
