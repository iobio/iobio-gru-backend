const process = require('process');
const spawn = require('child_process').spawn;
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');

const app = new Koa();
const router = new Router();

const aliases = {
  'samtools': 'anderspitman/samtools',
};

function run(command, args) {
  const child = spawn('docker', ['run', '--rm', aliases[command]].concat(args));
  return child;
}

router.get('/getAlignmentHeader', (ctx, next) => {
  const url = decodeURIComponent(ctx.query.url);

  const child = run('samtools', ['view', '-H', url]);

  ctx.body = child.stdout;

  console.log(JSON.stringify(ctx.query, null, 2));
});

router.get('/getAlignment', (ctx, next) => {
  const url = decodeURIComponent(ctx.query.url);

  const child = run('samtools', ['view', url, '18']);

  child.stderr.pipe(process.stderr);

  ctx.body = child.stdout;
});

app
  .use(cors())
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(9001);
