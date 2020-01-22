const spawn = require('child_process').spawn;

async function mktemp() {
  const proc = spawn('mktemp');

  let data = '';
  proc.stdout.on('data', (chunk) => {
    data += chunk;
  });

  return new Promise((resolve, reject) => {
    proc.on('exit', (e) => {
      if (e === 0) {
        resolve(data.slice(0, data.length - 1));
      }
      else {
        reject(e);
      }
    });
  });
}

//mktemp().then((tmp) => {
//  console.log(tmp, tmp.length);
//}).catch((e) => {
//  console.error(e);
//});

module.exports = {
  mktemp
};
