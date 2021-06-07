const path = require('path');

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

module.exports = {
  parseArgs,
  dataPath,
};
