import path from 'node:path';

function semverLess(a, b) {
  const pa = a.split('.').map(Number);
  const pb = b.split('.').map(Number);

  for (let i = 0; i < 3; i++) {
    if ((pa[i] || 0) !== (pb[i] || 0)) {
      return (pa[i] || 0) < (pb[i] || 0);
    }
  }
  return false;
}

function genRegionStr(region) {
  return region.refName + ':' + region.start + '-' + region.end;
}

function escapeRegExp(string) {
  return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); // $& means the whole matched string
}

function replaceAll(str, find, replace) {
  return str.replace(new RegExp(escapeRegExp(find), 'g'), replace);
}

function parseArgs() {
  const args = process.argv
    .slice(2)
    .map(arg => arg.split('='))
    .reduce((args, [value, key]) => {
        args[value] = key;
        return args;
    }, {});

  return args;
}

const args = parseArgs();

let dataDir = './data';
if (args['--data-dir']) {
  dataDir = args['--data-dir'];
}

function dataPath(name) {
  const absPath = path.resolve(path.join(dataDir, name));
  return absPath;
}

export {
  semverLess,
  genRegionStr,
  parseArgs,
  dataPath,
  replaceAll,
};
