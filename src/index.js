const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const path = require('path');
const { run } = require('./process.js');
const spawn = require('child_process').spawn;
const process = require('process');


const app = new Koa();
const router = new Router();


router.get('/', async (ctx) => {
  ctx.body = "<h1>I be healthful</h1>";
});

router.get('/alignmentHeader', async (ctx) => {
  await handle(ctx, 'alignmentHeader.sh', [ctx.query.url]);
});

router.get('/baiReadDepth', async (ctx) => {
  await handle(ctx, 'baiReadDepth.sh', [ctx.query.url]);
});

router.get('/craiReadDepth', async (ctx) => {
  await handle(ctx, 'craiReadDepth.sh', [ctx.query.url]);
});

router.get('/alignmentStatsStream', async (ctx) => {

  const regions = JSON.parse(ctx.query.regions);

  const samtoolsRegions = regions.map(function(d) { return d.name+ ":"+ d.start + '-' + d.end;}).join(' ');
  const bamstatsRegions = JSON.stringify(regions.map(function(d) { return {start:d.start,end:d.end,chr:d.name};}));

  await handle(ctx, 'alignmentStatsStream.sh', [ctx.query.url, samtoolsRegions, ctx.query.indexUrl, bamstatsRegions]);
});


async function handle(ctx, scriptName, args) {
  try {
    const scriptPath = path.join(__dirname, '../scripts', scriptName);
    const proc = await run(scriptPath, args);
    ctx.body = proc.stdout;
  }
  catch (e) {
    console.error(e);
    ctx.status = 500;
    ctx.body = e;
  }
}

app
  .use(cors())
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(9001);
