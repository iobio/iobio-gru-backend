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

// bam.iobio endpoints
//
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


// gene.iobio endpoints
router.get('/variantHeader', async (ctx) => {
  await handle(ctx, 'variantHeader.sh', [ctx.query.url]);
});

router.get('/vcfReadDepth', async (ctx) => {
  await handle(ctx, 'vcfReadDepth.sh', [ctx.query.url]);
});

router.get('/alignmentCoverage', async (ctx) => {

  const url = ctx.query.url;
  const indexUrl = ctx.query.indexUrl;
  const samtoolsRegion = JSON.parse(ctx.query.samtoolsRegion);
  const maxPoints = ctx.query.maxPoints;
  const coverageRegions = JSON.parse(ctx.query.coverageRegions);

  const samtoolsRegionArg = samtoolsRegion.refName + ':' + samtoolsRegion.start + '-' + samtoolsRegion.end;
  const spanningRegionArg = "-r " + samtoolsRegion.refName + ':' + samtoolsRegion.start + ':' + samtoolsRegion.end;

  const coverageRegionsArg = coverageRegions.length === 0 ? '' :
    "-p " + coverageRegions
      .filter(d => d.name && d.start && d.end)
      .map(d => d.name + ":" + d.start + ':' + d.end)
      .join(' ');

  const maxPointsArg = "-m " + maxPoints;

  const args = [url, indexUrl, samtoolsRegionArg, maxPointsArg, spanningRegionArg, coverageRegionsArg];

  await handle(ctx, 'alignmentCoverage.sh', args, { ignoreStderr: true });
});

async function handle(ctx, scriptName, args, options) {
  try {
    const scriptPath = path.join(__dirname, '../scripts', scriptName);
    const proc = await run(scriptPath, args, options ? options : {});
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
