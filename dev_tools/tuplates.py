# tuplate_start(https://cdn.jsdelivr.net/gh/anderspitman/tuplates@v0.3.0/tuplates.py)

import os, json
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

    with open ('tuplates_config.json') as f:
        config = parse_config(json.load(f))

    for (dirpath, dirnames, filenames) in os.walk("./"):

        del_list = [ (i, x) for (i, x) in enumerate(dirnames) if x in config['exclude'] ]
        for item in del_list:
            del dirnames[item[0]]

        for filename in filenames:
            path = os.path.join(dirpath, filename)

            _, ext = os.path.splitext(path)

            if config['extensions'] and ext.lower() not in config['extensions']:
                continue

            process_file(path)

# tuplate_end
