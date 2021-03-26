# tuplate_star(https://cdn.jsdelivr.net/gh/anderspitman/tuplates@v0.5.0/tuplates.py)

import os, json, argparse
from urllib import request

# These need to be built with concatenation so that tuplates.py can be run
# on itself (ie to update). Otherwise it falsely detects tuplates in the
# parsing code below.
start_str = 'tuplate_' + 'start('
end_str = 'tuplate_' + 'end'

def process_file(path):

    print ("Process", path)

    modified = False
    replacing = False

    out = ''

    with open(path, 'r') as f:
        try:
            for line in f:
                if start_str in line:
                    if replacing:
                        raise Exception("Unexpected " + start_str)
                    replacing = True

                    location = line.split(start_str)[1].split(')')[0]
                    tuplate = get_tuplate(location)

                    out += line
                    out += tuplate
                    modified = True

                elif end_str in line:
                    if not replacing:
                        raise Exception("Unexpected " + end_str)
                    replacing = False

                    out += line
                elif replacing:
                    pass
                else:
                    out += line
        except:
            print("Failed to process {}. Is it a text file?".format(path))
            return

    # reopen file to overwrite it
    if modified:
        with open(path, 'w') as f:
            f.write(out)

def get_tuplate(location):
    print ("    Get tuplate", location)
    if location.startswith('http'):
        res = request.urlopen(location)
        return res.read().decode('utf-8')
    else:
        with open(location) as f:
            return f.read()

def parse_config(config):
    keys = ['exclude', 'extensions']
    for key in keys:
        if key not in config:
            config[key] = []
    return config

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--config', default='tuplates_config.json', help='Config path')
    parser.add_argument('-t','--file-type', action='append', help='Add file type')
    parser.add_argument('-e','--exclude', action='append', help='Exclude')
    args = parser.parse_args()

    try:
        with open (args.config) as f:
            config = parse_config(json.load(f))
    except:
        config = parse_config({})

    if args.file_type:
        for ft in args.file_type:
            if ft not in config['extensions']:
                config['extensions'].append(ft)

    if args.exclude:
        for e in args.exclude:
            if e not in config['exclude']:
                config['exclude'].append(e)

    for (dirpath, dirnames, filenames) in os.walk("./"):

        del_list = [ (i, x) for (i, x) in enumerate(dirnames) if x in config['exclude'] ]

        # Go in reverse order since we need to delete in place and don't want
        # elements to shift as things are deleted.
        for item in reversed(del_list):
            del dirnames[item[0]]

        for filename in filenames:
            path = os.path.join(dirpath, filename)

            _, ext = os.path.splitext(path)
            # skip period
            ext = ext[1:]

            if config['extensions'] and ext.lower() not in config['extensions']:
                continue

            process_file(path)

# tuplate_en
