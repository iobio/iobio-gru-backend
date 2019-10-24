const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const logger = require('koa-logger');
const bodyParser = require('koa-bodyparser');
const path = require('path');
const { run } = require('./process.js');
const spawn = require('child_process').spawn;
const process = require('process');

let port = 9001;
if (process.argv[2]) {
  port = process.argv[2];
}

const app = new Koa();
const router = new Router();


const dataDir = './data';
function dataPath(name) {
  const absPath = path.resolve(path.join(dataDir, name));
  return absPath;
}

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

  const params = JSON.parse(ctx.request.body);
  console.log(JSON.stringify(params, null, 2));

  const url = params.url;
  const indexUrl = params.indexUrl;
  const samtoolsRegion = params.samtoolsRegion;
  const maxPoints = params.maxPoints;
  const coverageRegions = params.coverageRegions;

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
  const refFastaFile = dataPath(decodeURIComponent(ctx.query.refFastaFile));

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

  const q = ctx.query;
  console.log(JSON.stringify(q, null, 2));

  const tbiUrl = q.tbiUrl ? q.tbiUrl : '';
  const contigStr = genContigFileStr(JSON.parse(q.refNames));
  const regionStr = genRegionsStr(JSON.parse(q.regions));
  const vcfSampleNamesStr = JSON.parse(q.vcfSampleNames).join("\n");
  const refFastaFile = dataPath(q.refFastaFile);
  const vepCacheDir = dataPath('vep-cache');
  const vepREVELFile = dataPath(q.vepREVELFile);
  const vepPluginDir = dataPath('vep-cache/Plugins');

  const gnomadUrl = q.gnomadUrl ? q.gnomadUrl : '';
  const gnomadRegionStr = q.gnomadRegionStr ? q.gnomadRegionStr : '';
  const gnomadHeaderFile = dataPath('gnomad_header.txt');

  const args = [
    q.vcfUrl, tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
    refFastaFile, q.genomeBuildName, vepCacheDir, vepREVELFile, q.vepAF,
    vepPluginDir, q.isRefSeq, q.hgvsNotation, q.getRsId, gnomadUrl,
    gnomadRegionStr, gnomadHeaderFile,
  ];

  //await handle(ctx, 'annotateVariants.sh', args);
  await handle(ctx, 'annotateVariants.sh', args, { ignoreStderr: true });
});

router.post('/freebayesJointCall', async (ctx) => {

  const params = JSON.parse(ctx.request.body);
  console.log(JSON.stringify(params, null, 2));

  const alignments = params.alignmentSources.map(aln => aln.bamUrl).join(',');
  const indices = params.alignmentSources.map(aln => aln.baiUrl).join(',');
  const region = genRegionStr(params.region);
  const vepREVELFile = dataPath(params.vepREVELFile);
  const refFastaFile = dataPath(params.refFastaFile);
  const contigStr = genContigFileStr(params.refNames);
  const samplesFileStr = params.sampleNames.join('\n');

  const fbArgs = params.fbArgs;
  const freebayesArgs = [];
  if (fbArgs) {
    for (var key in fbArgs) {
      var theArg = fbArgs[key];
      if (theArg.hasOwnProperty('argName')) {
        if (theArg.hasOwnProperty('isFlag') && theArg.isFlag == true) {
          if (theArg.value && theArg.value == true) {
              freebayesArgs.push(theArg.argName);
          }
        } else {
          if (theArg.value && theArg.value != '') {
            freebayesArgs.push(theArg.argName);
            freebayesArgs.push(theArg.value);
          }
        }

      }
    }
  }

  const useSuggestedVariants = params.fbArgs.useSuggestedVariants.value ? 'true' : '';

  console.log(freebayesArgs);
  const extraArgs = freebayesArgs;

  const args = [
    alignments, indices, region, refFastaFile, useSuggestedVariants,
    params.clinvarUrl, params.genomeBuildName, vepREVELFile, params.vepAF,
    params.isRefSeq, samplesFileStr, extraArgs,
  ];

  await handle(ctx, 'freebayesJointCall.sh', args, { ignoreStderr: true });
});

router.post('/clinvarCountsForGene', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  console.log(JSON.stringify(params, null, 2));

  const region = genRegionStr(params.region);
  const regions = params.regions;

  let regionParts = "";
  if (regions) {
    regions.forEach(function(region) {
      if (regionParts.length > 0) {
        regionParts += ",";
      }
      regionParts += region.start + "-" + region.end;
    })
  }

  const binLength = params.binLength ? params.binLength : '';

  const args = [
    params.clinvarUrl, region, binLength, regionParts,
  ];

  await handle(ctx, 'clinvarCountsForGene.sh', args);
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

function genRegionStr(region) {
  return region.refName + ':' + region.start + '-' + region.end;
}

function genRegionsStr(regions) {
  let regionStr = "";
  for (const region of regions) {
    regionStr += region.name + ':' + region.start + '-' + region.end + " ";
  }
  return regionStr;
}

app
  .use(logger())
  .use(cors({
    maxAge: 86400,
  }))
  .use(bodyParser({
    enableTypes: ['json', 'text'],
  }))
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(port);
