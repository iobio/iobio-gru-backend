const process = require('process');
const spawn = require('child_process').spawn;
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');


const COMMANDS = {
  'samtools': {
    dockerImage: 'anderspitman/samtools',
  },
};

function run(command, args) {
  const child = spawn('docker', ['run', '--rm', COMMANDS[command].dockerImage].concat(args));
  return child;
}

class RPCServer {
  constructor(config) {

    this._config = config;
    this._app = new Koa();
    this._router = new Router();

    for (const methodName in config) {
      this._router.get('/' + methodName, (ctx, next) => {

        console.log(JSON.stringify(ctx.query, null, 2));

        // TODO: error handling for invalid params
        const params = this._parseParams(config[methodName], ctx.query);
        console.log(params);
        const pipeline = config[methodName].pipeline(params);
        const command = pipeline[0];
        const child = run(command[0], command.slice(1));

        if (config[methodName].returnStream === false) {
          // TODO: do something different here?
          ctx.body = child.stdout;
        }
        else {
          ctx.body = child.stdout;
        }
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

const rpcConfig = {

  getAlignmentHeader: {
    params: {
      url: 'URL',
    },
    pipeline: (params) => {
      return [
        ['samtools', 'view', '-H', params.url],
      ];
    },
    returnStream: false,
  },

  getAlignment: {
    params: {
      url: 'URL',
      chr: 'String',
    },
    pipeline: (params) => {
      return [
        ['samtools', 'view', params.url, params.chr],
      ];
    },
  },

};


const server = new RPCServer(rpcConfig);
server.start(9001);
