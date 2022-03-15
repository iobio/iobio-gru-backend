const spawn = require('child_process').spawn;


async function run(path, args, options) {

  const proc = spawn(path, args);

  if (options && options.inStream) {
    options.inStream.pipe(proc.stdin);
  }

  let settled = false;

  return new Promise((resolve, reject) => {

    proc.stderr.setEncoding('utf8');

    function onStdout(data) {
      proc.stdout.removeListener('data', onStdout);
      proc.stdout.pause();
      proc.stdout.unshift(data);
      resolve(proc);
    }
    proc.stdout.on('data', onStdout);

    function onStderr(data) {
      console.error(data);
      if (options && !options.ignoreStderr) {
        if (!settled) {
          settled = true;
          reject(new Error(data));
        }
      }
    }
    proc.stderr.on('data', onStderr);

    function onExit(e) {

      if (!settled) {
        settled = true;

        if (e === 0) {
          resolve(proc);
        }
        else {
          reject(new Error("return code: " + e));
        }
      }
    }
    proc.on('exit', onExit);

  });
}


module.exports = {
  run,
};
