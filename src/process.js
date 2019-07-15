const spawn = require('child_process').spawn;


// The purpose of this wrapper is to detect error conditions on the child
// process, before yielding the process to be piped. The basic idea is to
// listen to stdout and stderr, and watch for the exit event. If stdout
// receives data before stderr and/or exit, then the process is provided to
// the caller, otherwise the promise is rejected.
async function run(path, args, options) {

  const proc = spawn(path, args);

  if (options && options.inStream) {
    options.inStream.pipe(proc.stdin);
  }

  let settled = false;

  return new Promise((resolve, reject) => {

    //proc.stdout.setEncoding('utf8');
    proc.stderr.setEncoding('utf8');

    //const BUF_SIZE = 1024;
    //let buf = "";

    function onStdout(data) {

      //buf += data;

      //if (buf.length > BUF_SIZE) {
        proc.stdout.removeListener('data', onStdout);
        proc.stdout.pause();
        proc.stdout.unshift(data);
        //proc.stdout.unshift(buf);
        //buf = "";
        resolve(proc);
      //}
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
          // Just in case the overall length is less than BUF_SIZE, make sure
          // it gets passed on.
          //if (buf.length > 0) {
          //  proc.stdout.pause();
          //  console.log("unshift");
          //  proc.stdout.unshift(buf);
          //}
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
