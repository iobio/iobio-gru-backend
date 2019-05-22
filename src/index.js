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

    this._app = new Koa();
    this._router = new Router();

    for (const methodName in config) {
      this._router.get('/' + methodName, (ctx, next) => {

        const pipeline = config[methodName].pipeline(ctx.query);
        const command = pipeline[0];
        const child = run(command[0], command.slice(1));
        console.log(JSON.stringify(ctx.query, null, 2));

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
    },
    pipeline: (params) => {
      return [
        ['samtools', 'view', params.url, '18'],
      ];
    },
  },
};


const server = new RPCServer(rpcConfig);
server.start(9001);
