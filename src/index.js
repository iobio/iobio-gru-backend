const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const logger = require('koa-logger');
const bodyParser = require('koa-bodyparser');
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
//
router.get('/variantHeader', async (ctx) => {
  await handle(ctx, 'variantHeader.sh', [ctx.query.url, ctx.query.indexUrl]);
});

router.get('/vcfReadDepth', async (ctx) => {
  await handle(ctx, 'vcfReadDepth.sh', [ctx.query.url]);
});

router.post('/alignmentCoverage', async (ctx) => {

  console.log(JSON.stringify(ctx.request.body, null, 2));

  const url = ctx.request.body.url;
  const indexUrl = ctx.request.body.indexUrl;
  const samtoolsRegion = ctx.request.body.samtoolsRegion;
  const maxPoints = ctx.request.body.maxPoints;
  const coverageRegions = ctx.request.body.coverageRegions;

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

router.get('/geneCoverage', async (ctx) => {

  const url = ctx.query.url;
  const indexUrl = ctx.query.indexUrl;
  const refName = ctx.query.refName;
  const geneName = ctx.query.geneName;
  const regionStart = ctx.query.regionStart;
  const regionEnd = ctx.query.regionEnd;
  const regions = JSON.parse(ctx.query.regions);

  let regionStr = "#" + geneName + "\n";
  regions.forEach(function(region) {
    regionStr += refName + ":" + region.start + "-" + region.end + "\n";
  });

  const samtoolsRegionArg = refName + ':' + regionStart + '-' + regionEnd;

  const args = [url, indexUrl, samtoolsRegionArg, regionStr];

  await handle(ctx, 'geneCoverage.sh', args);
});

router.get('/normalizeVariants', async (ctx) => {

  const vcfUrl = ctx.query.vcfUrl;
  const tbiUrl = ctx.query.tbiUrl;
  const refName = ctx.query.refName;
  const regions = JSON.parse(ctx.query.regions);
  const contigStr = decodeURIComponent(ctx.query.contigStr);
  const refFastaFile = './data/' + decodeURIComponent(ctx.query.refFastaFile);

  let regionParm = "";
  regions.forEach(function(region) {
    if (regionParm.length > 0) {
      regionParm += " ";
    }
    regionParm += region.refName + ":" + region.start + "-" + region.end;
  });

  const args = [vcfUrl, tbiUrl, refName, regionParm, contigStr, refFastaFile];

  await handle(ctx, 'normalizeVariants.sh', args);
});

router.get('/annotateVariants', async (ctx) => {

  console.log(ctx.query);
  const q = ctx.query;

  const contigStr = genContigFileStr(JSON.parse(q.refNames));
  const regionStr = genRegionStr(JSON.parse(q.regions));
  const vcfSampleNamesStr = JSON.parse(q.vcfSampleNames).join("\n");
  const vepREVELFile = './data/' + q.vepREVELFile;
  const refFastaFile = './data/' + q.refFastaFile;

  const args = [
    q.vcfUrl, q.tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
    refFastaFile, q.genomeBuildName, vepREVELFile, q.vepAF, q.isRefSeq,
    q.hgvsNotation, q.getRsId,
  ];

  console.log(args);
  //await handle(ctx, 'annotateVariants.sh', args);
  await handle(ctx, 'annotateVariants.sh', args, { ignoreStderr: true });
});






// genepanel endpoints
//
router.get('/clinphen', async (ctx) => {

  const args = [ctx.query.notes];

  await handle(ctx, 'clinphen.sh', args);
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

function genContigFileStr(refNames) {
  let contigStr = "";
  for (const ref of refNames) {
    contigStr += "##contig=<ID=" + ref + ">\n";
  }
  return contigStr;
}

function genRegionStr(regions) {
  let regionStr = "";
  for (const region of regions) {
    regionStr += region.name + ':' + region.start + '-' + region.end + " ";
  }
  return regionStr;
}

app
  .use(logger())
  .use(cors())
  .use(bodyParser())
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(9001);
