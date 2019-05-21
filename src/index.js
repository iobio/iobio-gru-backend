const process = require('process');
const spawn = require('child_process').spawn;
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');

const app = new Koa();
const router = new Router();

router.get('/getAlignmentHeader', (ctx, next) => {
  const alignmentUrl = decodeURIComponent(ctx.query.alignmentUrl);

  const child = spawn('samtools', ['view', '-H', alignmentUrl]);

  ctx.body = child.stdout;

  console.log(JSON.stringify(ctx.query, null, 2));
});

router.get('/getAlignment', (ctx, next) => {
  const alignmentUrl = decodeURIComponent(ctx.query.alignmentUrl);

  const child = spawn('samtools', ['view', alignmentUrl, '18']);

  child.stderr.pipe(process.stderr);

  ctx.body = child.stdout;
});

app
  .use(cors())
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(9001);
