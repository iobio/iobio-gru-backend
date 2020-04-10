const path = require('path');

const dataDir = './data';
function dataPath(name) {
  const absPath = path.resolve(path.join(dataDir, name));
  return absPath;
}

module.exports = {
  dataPath,
};
