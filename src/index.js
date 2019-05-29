const spawn = require('child_process').spawn;
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');


function run(tool, args) {
  const toolPath = './tools/' + tool;
  const child = spawn(toolPath, args);
  return child;
}

class RPCServer {
  constructor(config) {

    this._config = config;
    this._app = new Koa();
    this._router = new Router();

    for (const methodName in config) {
      this._router.get('/' + methodName, (ctx, next) => {

        // TODO: error handling for invalid params
        const params = this._parseParams(config[methodName], ctx.query);
        const pipeline = config[methodName].pipeline(params);

        //const command = pipeline[0];

        let prev;
        let child;

        for (const stage of pipeline) {
          child = run(stage[0], stage.slice(1));

          if (prev) {
            prev.stdout.pipe(child.stdin);
          }

          prev = child;
        }

        ctx.body = child.stdout;
      });
    }

  }

  start(port) {
    this._app
      .use(cors())
      .use(this._router.routes())
      .use(this._router.allowedMethods())
      .listen(port);
  }

  _parseParams(method, queryParams) {
    const params = {};

    for (const key in method.params) {
      const type = method.params[key];

      switch (type) {
        case 'URL':
          params[key] = decodeURIComponent(queryParams[key]);
          break;
        case 'String':
          params[key] = queryParams[key];
          break;
        default:
          throw new Error("Invalid type: " + type);
          break;
      }
    }

    return params;
  }
}

module.exports = {
  RPCServer,
};
