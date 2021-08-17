#!/usr/bin/env python3

# tuplate_start(https://cdn.jsdelivr.net/gh/gemdrive/gemdrive-cli-py@v0.7.0/client.py)
import os, threading, queue, json, math, shutil, stat
from urllib import request, parse
from datetime import datetime


# Ensures default values are properly set if missing
def clean_gem_data(data):

    if 'size' not in data:
        data['size'] = 0

    return data


class GemDriveClient():

    def __init__(self, **kwargs):

        self.options = kwargs
        self.job_queue = queue.Queue(maxsize=8)

        for w in range(kwargs['num_workers']):
            threading.Thread(target=self.downloader, daemon=True).start()

    def sync(self, gemdrive_url, fs_dir):
        self.traverse(gemdrive_url, fs_dir, None)
        self.job_queue.join()

    def traverse(self, url, parent_dir, gem_data_in):

        depth = self.options['depth']
        token = self.options['token']

        gem_data = gem_data_in

        if gem_data is None:
            u = parse.urlparse(url)
            p = parse.quote(u.path)
            gem_url = u.scheme + '://' + u.netloc + '/gemdrive/index' + p + 'tree.json?depth=' + str(depth)

            if token is not None:
                gem_url += '&access_token=' + token

            try:
                res = request.urlopen(gem_url, timeout=5)
                body = res.read()
                gem_data = json.loads(body)
            except request.URLError as e:
                print("Timed out retrieving " + gem_url)
                raise

        gem_data = clean_gem_data(gem_data)

        if not os.path.isdir(parent_dir):
            print("Create", parent_dir)
            if self.options['dry_run']:
                # Early return for dry run because attempts to compare child
                # directories which don't exist will cause exceptions.
                return
            else:
                os.mkdir(parent_dir)

        if 'children' not in gem_data:
            return

        for child_name in gem_data['children']:
            child = gem_data['children'][child_name]
            child = clean_gem_data(child)
            child_url = url + child_name
            child_path = os.path.join(parent_dir, child_name)
            is_dir = child_url.endswith('/')
            if is_dir:
                self.traverse(child_url, child_path, child)
            else:
                self.job_queue.put((child_url, parent_dir, child))

        if self.options['delete']:
            with os.scandir(parent_dir) as it:
                for entry in it:
                    name = entry.name
                    if entry.is_dir():
                        name += '/'

                    if name not in gem_data['children']:
                        item_path = os.path.join(parent_dir, name)
                        print("Delete", item_path)

                        if not self.options['dry_run']:
                            if entry.is_dir():
                                shutil.rmtree(item_path)
                            else:
                                os.remove(item_path)

    def downloader(self):
        while True:
            url, parent_dir, gem_data = self.job_queue.get()

            self.handle_file(url, parent_dir, gem_data)

            self.job_queue.task_done()


    def handle_file(self, url, parent_dir, gem_data):

        token = self.options['token']

        name = os.path.basename(url)
        path = os.path.join(parent_dir, name)

        try:
            stats = os.stat(path)
        except:
            stats = None

        size = 0
        mod_time = ''
        dest_is_exe = False
        if stats:
            size = stats.st_size
            mod_time = datetime.utcfromtimestamp(stats.st_mtime).replace(microsecond=0).isoformat() + 'Z'
            dest_is_exe = stats.st_mode & 0o111 != 0


        utc_dt = datetime.strptime(gem_data['modTime'], '%Y-%m-%dT%H:%M:%SZ')
        mtime = math.floor((utc_dt - datetime(1970, 1, 1)).total_seconds())

        needs_update = size != gem_data['size'] or mod_time != gem_data['modTime']

        src_is_exe = 'isExecutable' in gem_data and gem_data['isExecutable']

        if src_is_exe != dest_is_exe:
            needs_update = True

        if needs_update:
            print("Sync", path)

            if not self.options['dry_run']:
                file_url = url

                if token is not None:
                    file_url += '?access_token=' + token

                u = parse.urlparse(file_url)
                req_url = u.scheme + '://' + u.netloc + parse.quote(u.path) + '?' + u.query

                try:
                    res = request.urlopen(req_url, timeout=5)
                except request.URLError as e:
                    print("Timed out retrieving " + req_url)
                    raise

                with open(path, 'wb') as f:
                    while True:
                        chunk = res.read(4096)
                        if not chunk:
                            break
                        f.write(chunk)
                stats = os.stat(path)

                if stats.st_size != gem_data['size']:
                    print("Sizes don't match", url)

                os.utime(path, (stats.st_atime, mtime))

                if src_is_exe and not dest_is_exe:
                    os.chmod(path, stats.st_mode | 0o111)
# tuplate_end


import os, argparse


def dir_name(path):
    return os.path.basename(path[:-1])

if __name__ == '__main__':

    cwd = os.getcwd()

    parser = argparse.ArgumentParser()
    parser.add_argument('--gru-version', help='gru version to download', default='1.2.1')
    parser.add_argument('--out-dir', help='Output directory', default=cwd)
    parser.add_argument('--dry-run', help='Enable dry run mode. No changes will be made to destination',
            default=False, action='store_true')
    args = parser.parse_args()

    client = GemDriveClient(depth=0, token=None, verbose=False, dry_run=args.dry_run, delete=True,
            num_workers=8)

    url = 'https://gemdrive.iobio.io/gru-' + args.gru_version + '/'

    if args.out_dir == cwd:
        args.out_dir = os.path.join(cwd, dir_name(url))

    client.sync(url, args.out_dir)
