const spawn = require('child_process').spawn;
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');

const { run } = require('./process');


class RPCServer {
  constructor(config) {

    this._config = config;
    this._app = new Koa();
    this._router = new Router();

    this._router.get('/call', async (ctx, next) => {

      const pipeline = JSON.parse(ctx.query.pipeline);
      console.log(pipeline);

      // There is an attempt at some error handling in the block below.
      // The basic algorithm is listen for stderr (can be ignored) and
      // exit with code !== 0 for all processes in the pipeline. Also listen
      // for stdout of the last process. If stderr or exit !== 0 occurs
      // before a specified amount of stdout buffer is received, then
      // consider it an error condition.
      // If an error occurs after we've started streaming, we're currently
      // not handling that. I think we'd either need to add HTTP trailers
      // or switch to WebSockets.
      const promise = new Promise((resolve, reject) => {

        let aborted = false;

        function onStderr(data) {
          if (!aborted) {
            aborted = true;

            reject(data);
          }
        }

        function onExit(exitCode, command) {
          if (!aborted) {
            aborted = true;

            // TODO: if already started streaming, either send an HTTP
            // trailer here or switch to WebSockets so we can report the
            // error.
            if (exitCode !== 0) {
              reject(`ERROR: Process exited with code ${exitCode}: ${command}`);
            }
          }
        }


        let prev;
        let child;

        for (const stage of pipeline) {
          console.log(new Date());

          const name = stage[0];
          const args = stage[1];
          const options = stage[2];

          child = spawn('./tools/' + name, args);

          if (options && options.ignoreStderr !== true) {
            child.stderr.setEncoding('utf8');
            child.stderr.on('data', onStderr);
          }

          child.on('exit', (exitCode) => {
            onExit(exitCode, JSON.stringify(stage));
          });

          if (prev) {
            prev.stdout.pipe(child.stdin);
          }

          prev = child;
        }

        const BUF_SIZE = 1024;
        let buf = "";
        function onStdout(data) {

          buf += data;

          if (buf.length > BUF_SIZE) {
            child.stdout.removeListener('data', onStdout);
            child.stdout.pause();
            child.stdout.unshift(buf);
            resolve(child.stdout);
          }
        }

        child.stdout.on('data', onStdout);
      });

      try {
        const stream = await promise;
        ctx.body = stream;
      }
      catch (e) {
        ctx.status = 500;
        ctx.body = e;
      }
    });


  }

  start(port) {
    this._app
      .use(cors())
      .use(this._router.routes())
      .use(this._router.allowedMethods())
      .listen(port);
  }
}


module.exports = {
  RPCServer,
};
