const spawn = require('child_process').spawn;


// The purpose of this wrapper is to detect error conditions on the child
// process, before yielding the process to be piped. The basic idea is to
// listen to stdout and stderr, and watch for the exit event. If stdout
// receives data before stderr and/or exit, then the process is provided to
// the caller, otherwise the promise is rejected.
async function run(path, args, options) {

  const proc = spawn(path, args);

  if (options && options.inStream) {
    console.log("pipeit");
    options.inStream.pipe(proc.stdin);
  }

  let settled = false;

  return new Promise((resolve, reject) => {

    //proc.stdout.setEncoding('utf8');
    //proc.stderr.setEncoding('utf8');

    function onStdout(data) {
      console.log("stdout");
      //console.log(data);
      settled = true;

      proc.stdout.removeListener('data', onStdout);
      proc.stderr.removeListener('data', onStderr);
      proc.removeListener('exit', onExit);

      // put the data back so the caller will read it. Need to pause since
      // we're currently streaming.
      proc.stdout.pause();
      proc.stdout.unshift(data);

      resolve(proc);
    }
    proc.stdout.on('data', onStdout);

    function onStderr(data) {
      //console.log("stderr");
      //console.log(data);
      if (options && !options.ignoreStderr) {
        console.log("stderr");
        if (!settled) {
          settled = true;
          reject(new Error(data));
        }
      }
    }
    proc.stderr.on('data', onStderr);

    function onExit(e) {
      console.log("exit");
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
