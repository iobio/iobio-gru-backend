const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const logger = require('koa-logger');
const bodyParser = require('koa-bodyparser');
const mount = require('koa-mount');
const serve = require('koa-static');
const path = require('path');
const { run } = require('./process.js');
const { mktemp } = require('./mktemp.js');
const spawn = require('child_process').spawn;
const process = require('process');
const gene2PhenoRouter = require('./gene2pheno.js');
const geneInfoRouter = require('./geneinfo.js');
const genomeBuildRouter = require('./genomebuild.js');
const hpoRouter = require('./hpo.js');
const { dataPath } = require('./utils.js');
const fs = require('fs');

let port = 9001;
if (process.argv[2]) {
  port = process.argv[2];
}

const app = new Koa();
const router = new Router();

router.use('/geneinfo', geneInfoRouter.routes(), geneInfoRouter.allowedMethods());
router.use('/gene2pheno', gene2PhenoRouter.routes(), gene2PhenoRouter.allowedMethods());
router.use('/genomebuild', genomeBuildRouter.routes(), genomeBuildRouter.allowedMethods());
router.use('/hpo', hpoRouter.routes(), hpoRouter.allowedMethods());

const staticServer = new Koa();
staticServer.use(serve(path.join(__dirname, '../static')));

router.get('/', async (ctx) => {
  ctx.body = "<h1>I be healthful</h1>";
});

router.post('/viewAlignments', async (ctx) => {
  const params = JSON.parse(ctx.request.body);

  const args = [params.url];

  if (params.regions) {
    const samtoolsRegions = genRegionsStr(params.regions);
    args.push(samtoolsRegions);
  }

  await handle(ctx, 'viewAlignments.sh', args);
});

// bam.iobio endpoints
//
// TODO: remove get in favor of post
router.get('/alignmentHeader', async (ctx) => {
  await handle(ctx, 'alignmentHeader.sh', [ctx.query.url]);
});
router.post('/alignmentHeader', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  await handle(ctx, 'alignmentHeader.sh', [params.url]);
});

router.get('/baiReadDepth', async (ctx) => {
  await handle(ctx, 'baiReadDepth.sh', [ctx.query.url]);
});
router.post('/baiReadDepth', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  await handle(ctx, 'baiReadDepth.sh', [params.url]);
});

router.get('/craiReadDepth', async (ctx) => {
  await handle(ctx, 'craiReadDepth.sh', [ctx.query.url]);
});
router.post('/craiReadDepth', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  await handle(ctx, 'craiReadDepth.sh', [params.url]);
});

router.get('/alignmentStatsStream', async (ctx) => {

  const regions = JSON.parse(ctx.query.regions);

  const samtoolsRegions = regions.map(function(d) { return d.name+ ":"+ d.start + '-' + d.end;}).join(' ');
  const bamstatsRegions = JSON.stringify(regions.map(function(d) { return {start:d.start,end:d.end,chr:d.name};}));

  await handle(ctx, 'alignmentStatsStream.sh', [ctx.query.url, samtoolsRegions, ctx.query.indexUrl, bamstatsRegions]);
});
router.post('/alignmentStatsStream', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

  const samtoolsRegions = genRegionsStr(params.regions);
  const bamstatsRegions = JSON.stringify(params.regions.map(function(d) { return {start:d.start,end:d.end,chr:d.name};}));

  await handle(ctx, 'alignmentStatsStream.sh', [params.url, samtoolsRegions, params.indexUrl, bamstatsRegions]);
});






// gene.iobio & oncogene & cohort-gene endpoints
//
// TODO: test post version and delete get above 
router.post('/variantHeader', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    await handle(ctx, 'variantHeader.sh', [params.url, params.indexUrl]);
});

router.post('/getChromosomes', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));
    await handle(ctx, 'getChromosomes.sh', [params.url, params.indexUrl]);
});

router.get('/vcfReadDepth', async (ctx) => {
  await handle(ctx, 'vcfReadDepth.sh', [ctx.query.url]);
});
// TODO: test post version and delete above
router.post('/vcfReadDepth', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    await handle(ctx, 'vcfReadDepth.sh', [params.url]);
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
      .join(',');
  console.log('REGION ARG GOING INTO ALIGNMENTCOVERAGE.SH');
  console.log(JSON.stringify(coverageRegionsArg));

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
// TODO: test post version and delete above
router.post('/geneCoverage', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

    // Pull from passed params
    const url = params.url;
    const indexUrl = params.indexUrl ? params.indexUrl : '';
    const refName = params.refName;
    const geneName = params.geneName;
    const regionStart = params.regionStart;
    const regionEnd = params.regionEnd;
    const regions = params.regions;

    // Format params
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
// TODO: test post version and delete above
router.post('/normalizeVariants', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

    // Pull from passed params
    const vcfUrl = params.vcfUrl;
    const tbiUrl = params.tbiUrl;
    const refName = params.refName;
    const regions = params.regions;
    const contigStr = decodeURIComponent(params.contigStr);
    const refFastaFile = dataPath(decodeURIComponent(params.refFastaFile));

    // Format params
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
    gnomadRegionStr, gnomadHeaderFile, q.decompose
  ];

  //await handle(ctx, 'annotateVariants.sh', args);
  await handle(ctx, 'annotateVariants.sh', args, { ignoreStderr: true });
});
router.post('/annotateVariants', async (ctx) => {

    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
    const contigStr = genContigFileStr(params.refNames);
    const regionStr = genRegionsStr(params.regions);
    const vcfSampleNamesStr = params.vcfSampleNames.join("\n");
    const refFastaFile = dataPath(params.refFastaFile);
    const vepCacheDir = dataPath('vep-cache');
    const vepREVELFile = dataPath(params.vepREVELFile);
    const vepPluginDir = dataPath('vep-cache/Plugins');

    const gnomadUrl = params.gnomadUrl ? params.gnomadUrl : '';
    const gnomadRegionStr = params.gnomadRegionStr ? params.gnomadRegionStr : '';
    const gnomadHeaderFile = dataPath('gnomad_header.txt');

    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
        refFastaFile, params.genomeBuildName, vepCacheDir, vepREVELFile, params.vepAF,
        vepPluginDir, params.isRefSeq, params.hgvsNotation, params.getRsId, gnomadUrl,
        gnomadRegionStr, gnomadHeaderFile, params.decompose
    ];

    await handle(ctx, 'annotateVariants.sh', args, { ignoreStderr: true });
});

router.post('/annotateEnrichmentCounts', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
    const contigStr = genContigFileStr(params.refNames);
    const regionStr = genRegionsStr(params.regions);
    const refFastaFile = dataPath(params.refFastaFile);
    const filterArgs = params.filterArgs ? params.filterArgs : '';
    const experStr = params.expIdString ? params.expIdString : '';
    const controlStr = params.controlIdString ? params.controlIdString : '';

    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr,
        refFastaFile, params.filterArgs,
        experStr, controlStr
    ];

    await handle(ctx, 'annotateEnrichmentCounts.sh', args, { ignoreStderr: true });
});

router.post('/getSomaticVariants', async (ctx) => {
  
  const params = JSON.parse(ctx.request.body);
  console.log(JSON.stringify(params, null, 2));

  const args = [params.vcfUrl, params.qualCutoff, params.totalReadCutoff, params.normalCountCutoff, params.tumorCountCutoff, params.normalAfCutoff, params.tumorAfCutoff, params.normalSampleIdx, params.totalSampleNum];
  
  await handle(ctx, 'getSomaticVariants.sh', args, { ignoreStderr: false });
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

  const vepCacheDir = dataPath('vep-cache');
  const vepPluginDir = dataPath('vep-cache/Plugins');


  const gnomadUrl = params.gnomadUrl ? params.gnomadUrl : '';
  const gnomadRegionStr = params.gnomadRegionStr ? params.gnomadRegionStr : '';
  const gnomadHeaderFile = dataPath('gnomad_header.txt');
  const decompose = params.decompose ? params.decompose : '';


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
    params.isRefSeq, samplesFileStr, extraArgs, vepCacheDir, vepPluginDir,
    gnomadUrl, gnomadRegionStr, gnomadHeaderFile, decompose
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
  const annotationMode = params.annotationMode ? params.annotationMode : '';
  const requiresVepService = params.requiresVepService ? params.requiresVepService : false;
  const vepArgs = params.vepArgs ? params.vepArgs : '';

  const args = [
    params.clinvarUrl, region, binLength, regionParts, annotationMode, requiresVepService, vepArgs
  ];

  await handle(ctx, 'clinvarCountsForGene.sh', args);
});



// genepanel endpoints
//
router.get('/clinphen', async (ctx) => {
 
  const args = [ctx.query.notes];

  await handle(ctx, 'clinphen.sh', args);
});

router.get('/phenotypeExtractor', async (ctx) => { 

  const args = [ctx.query.notes];

  await handle(ctx, 'phenotypeExtractor.sh', args);
});

router.post('/clinReport', async (ctx) => { 	 
  // Copy the data into a temporary file and then pass the path. It was failing
  // before, I'm pretty sure because the file was too large (~3MB) to pass
  // through the child spawing interface.
  const tmpFilePath = await mktemp();
  await fs.promises.writeFile(tmpFilePath, ctx.request.body);
  const args = [tmpFilePath];
  await handle(ctx, 'clinReport.sh', args); 	
});



// oncogene endpoints
//
router.post('/getIdColumns', async (ctx) => {

    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const regionStr = genRegionsStr(params.regions);
    const args = [
        params.vcfUrl, regionStr
    ];

    await handle(ctx, 'getIdColumns.sh', args, { ignoreStderr: true });
});

router.post('/checkBamBai', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const args = [ params.url, params.indexUrl, params.region ];
    await handle(ctx, 'checkBamBai.sh', args, { ignoreStderr: true });
});

// vcf.iobio endpoints
router.post('/vcfStatsStream', async (ctx) => {

  const params = JSON.parse(ctx.request.body);
  console.log(params);

  const regionStr = genRegionsStr(params.regions);
  const contigStr = genContigFileStr(params.refNames);

  let sampleNamesStr = "";
  if (params.sampleNames) {
    sampleNamesStr = params.sampleNames.join('\n');
  }

  const args = [
    params.url, params.indexUrl, regionStr, contigStr, sampleNamesStr
  ];
  console.log(args);

  await handle(ctx, 'vcfStatsStream.sh', args, { ignoreStderr: true });
}); 


async function handle(ctx, scriptName, args, options) {
  try {
    const scriptPath = path.join(__dirname, '../scripts', scriptName);
    const proc = await run(scriptPath, args, options ? options : {});
    ctx.body = proc.stdout;
  }
  catch (e) {
    console.error(e);
    ctx.status = 400;
    ctx.body = e.toString();
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

    regionStr += region.name;

    if (region.start) {
      regionStr += ':' + region.start;

      if (region.end) {
        regionStr += '-' + region.end + " ";
      }
    }
  }
  return regionStr;
}

app
  .use(mount('/static', staticServer))
  .use(logger())
  .use(cors({
    origin: '*',
    maxAge: 86400,
  }))
  .use(bodyParser({
    enableTypes: ['json', 'text'],
    jsonLimit: '10mb',
    textLimit: '10mb',
  }))
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(port);
