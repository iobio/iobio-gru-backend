const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const logger = require('koa-logger');
const bodyParser = require('koa-bodyparser');
const path = require('path');
const { run } = require('./process.js');
const { mktemp } = require('./mktemp.js');
const spawn = require('child_process').spawn;
const process = require('process');
const gene2PhenoRouter = require('./gene2pheno.js');
const geneInfoRouter = require('./geneinfo.js');
const genomeBuildRouter = require('./genomebuild.js');
const hpoRouter = require('./hpo.js');
const { parseArgs, dataPath } = require('./utils.js');
const fs = require('fs');
const { serveStatic } = require('./static.js');


const router = new Router();

router.use('/geneinfo', geneInfoRouter.routes(), geneInfoRouter.allowedMethods());
router.use('/gene2pheno', gene2PhenoRouter.routes(), gene2PhenoRouter.allowedMethods());
router.use('/genomebuild', genomeBuildRouter.routes(), genomeBuildRouter.allowedMethods());
router.use('/hpo', hpoRouter.routes(), hpoRouter.allowedMethods());

router.get('/static/*', async (ctx) => {

  // TODO: Hack to remove prefix when hosting app within gru because I don't
  // really understand how to do this properly with koa-router
  const reqPath = ctx.path.startsWith('/gru') ? ctx.path.slice(4) : ctx.path;

  const fsPath = path.join(__dirname, '..', reqPath);
  await serveStatic(ctx, fsPath);
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
router.post('/alignmentHeader', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  console.log('url going into alignment header: ' + params.url);
  await handle(ctx, 'alignmentHeader.sh', [params.url]);
});

router.post('/baiReadDepth', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  await handle(ctx, 'baiReadDepth.sh', [params.url]);
});

router.post('/craiReadDepth', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  await handle(ctx, 'craiReadDepth.sh', [params.url]);
});

router.post('/alignmentStatsStream', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

  const samtoolsRegions = genRegionsStr(params.regions);
  const bamstatsRegions = JSON.stringify(params.regions.map(function(d) { return {start:d.start,end:d.end,chr:d.name};}));

  const indexUrl = params.indexUrl ? params.indexUrl : '';

  await handle(ctx, 'alignmentStatsStream.sh', [
    params.url,
    indexUrl,
    samtoolsRegions,
    bamstatsRegions,
    dataPath(''),
  ], { ignoreStderr: true });
});


// gene.iobio & oncogene & cohort-gene endpoints
//
router.post('/variantHeader', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    const indexUrl = params.indexUrl ? params.indexUrl : '';
    await handle(ctx, 'variantHeader.sh', [params.url, indexUrl]);
});

router.post('/getChromosomes', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));
    const indexUrl = params.indexUrl ? params.indexUrl : '';
    await handle(ctx, 'getChromosomes.sh', [params.url, indexUrl]);
});

router.post('/vcfReadDepth', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    await handle(ctx, 'vcfReadDepth.sh', [params.url]);
});

router.post('/alignmentCoverage', async (ctx) => {

  const params = JSON.parse(ctx.request.body);
   console.log(JSON.stringify(params, null, 2));

  const url = params.url;
  const indexUrl = params.indexUrl ? params.indexUrl : '';
  const samtoolsRegion = params.samtoolsRegion;
  const maxPoints = params.maxPoints;
  const coverageRegions = params.coverageRegions;
  const qualityCutoff = params.qualityCutoff;
  const samtoolsRegionArg = samtoolsRegion.refName + ':' + samtoolsRegion.start + '-' + samtoolsRegion.end;
  const spanningRegionArg = samtoolsRegion.refName + ':' + samtoolsRegion.start + ':' + samtoolsRegion.end;

  const coverageRegionsArg = coverageRegions.length === 0 ? '' :
    coverageRegions
      .filter(d => d.name && d.start && d.end)
      .map(d => d.name + ":" + d.start + ':' + d.end)
      .join(',');

  const maxPointsArg = maxPoints;

  const args = [
    url, indexUrl, samtoolsRegionArg, maxPointsArg, spanningRegionArg,
    coverageRegionsArg, qualityCutoff, dataPath(''),
  ];

  await handle(ctx, 'alignmentCoverage.sh', args, { ignoreStderr: true });
});

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

    const dataDir = dataPath('');

    // Format params
    let regionStr = "#" + geneName + "\n";
    regions.forEach(function(region) {
        regionStr += refName + ":" + region.start + "-" + region.end + "\n";
    });
    const samtoolsRegionArg = refName + ':' + regionStart + '-' + regionEnd;
    const args = [url, indexUrl, samtoolsRegionArg, regionStr, dataDir];

    await handle(ctx, 'geneCoverage.sh', args, { ignoreStderr: true });
});

router.post('/normalizeVariants', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

    // Pull from passed params
    const vcfUrl = params.vcfUrl;
    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
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

router.post('/getClinvarVariants', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
    const contigStr = genContigFileStr(params.refNames);
    const regionStr = genRegionsStr(params.regions);
    const refFastaFile = dataPath(params.refFastaFile);
    const gnomadUrl = params.gnomadUrl ? params.gnomadUrl : '';
    const gnomadRegionStr = params.gnomadRegionStr ? params.gnomadRegionStr : '';
    const gnomadHeaderFile = dataPath('gnomad_header.txt');
    const gnomadRenameChr = params.gnomadRenameChr ? params.gnomadRenameChr : '';

    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr, refFastaFile, 
        params.genomeBuildName, gnomadUrl, gnomadRegionStr, 
        gnomadHeaderFile, gnomadRenameChr, params.clinSigFilterPhrase
    ];

    console.log(args);
   
    await handle(ctx, 'getClinvarVariants.sh', args, { ignoreStderr: true });
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
    const gnomadRenameChr = params.gnomadRenameChr ? params.gnomadRenameChr : '';


    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
        refFastaFile, params.genomeBuildName, vepCacheDir, vepREVELFile, params.vepAF,
        vepPluginDir, params.hgvsNotation, params.getRsId, gnomadUrl,
        gnomadRegionStr, gnomadHeaderFile, params.decompose, gnomadRenameChr

    ];

    console.log(args);
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

router.post('/annotateSomaticVariants', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  const vepCacheDir = dataPath('vep-cache');

  let refFastaFile = dataPath('references/GRCh37/human_g1k_v37_decoy_phix.fasta');
  if (params.genomeBuildName === 'GRCh38') {
    refFastaFile = dataPath('references/GRCh38/human_g1k_v38_decoy_phix.fasta');
  }
  
  const args = [params.vcfUrl, params.selectedSamplesStr, params.geneRegionsStr, params.somaticFilterPhrase, params.genomeBuildName, vepCacheDir, refFastaFile];
  
  await handle(ctx, 'annotateSomaticVariants.sh', args, { ignoreStderr: true });
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
    samplesFileStr, extraArgs, vepCacheDir, vepPluginDir,
    gnomadUrl, gnomadRegionStr, gnomadHeaderFile, decompose, dataPath(''),
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

    const regionStr = genRegionsStr(params.regions, ",");
    const args = [
        params.vcfUrl, regionStr
    ];

    console.log(regionStr);
    await handle(ctx, 'getIdColumns.sh', args, { ignoreStderr: true });
});

router.post('/checkBamBai', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    console.log(JSON.stringify(params, null, 2));

    const indexUrl = params.indexUrl ? params.indexUrl : '';

    const args = [ params.url, indexUrl, params.region, dataPath('') ];
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

  const indexUrl = params.indexUrl ? params.indexUrl : '';

  const args = [
    params.url, indexUrl, regionStr, contigStr, sampleNamesStr
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

function genRegionsStr(regions, delim = " ") {
  let regionStr = "";
  for (const region of regions) {

    regionStr += region.name;

    if (region.start) {
      regionStr += ':' + region.start;

      if (region.end) {
        regionStr += '-' + region.end + delim;
      }
    }
  }
  regionStr = regionStr.substring(0, regionStr.length - 1);
  return regionStr;
}




const args = parseArgs();

// This gives singularity images access to the data directory
process.env.SINGULARITY_BIND = dataPath('');

let toolDir = path.join(__dirname, '..', 'tool_bin');
if (args['--tools-dir']) {
  toolDir = args['--tools-dir'];
}
process.env.PATH = toolDir + ':' + process.env.PATH;

let port = 9001;
if (args['--port']) {
  port = Number(args['--port']);
}


let rootRouter = router;

// Allows a frontend app to be hosted in the same process. This is particularly
// useful if you want to run an entire app frontend and backend in a single
// docker container.
if (args['--app-dir']) {
  rootRouter = new Router();
  rootRouter.use('/gru', router.routes(), router.allowedMethods());

  rootRouter.get('/gru', async (ctx) => {
    ctx.body = "<h1>I be healthful</h1>";
  });

  rootRouter.get('/*', async (ctx, next) => {

    const fsPath = path.join(args['--app-dir'], ctx.path);
    await serveStatic(ctx, fsPath);

    // Bit of a hack. If previous attempt to serve didn't find the file,
    // default to the root index.html.
    if (ctx.status === 404) {
      const fsPath = path.join(args['--app-dir'], 'index.html');
      await serveStatic(ctx, fsPath);
    }
  });
}
else {
  rootRouter.get('/', async (ctx) => {
    ctx.body = "<h1>I be healthful</h1>";
  });
}

const server = new Koa();
server
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
  .use(rootRouter.routes())
  .use(rootRouter.allowedMethods())
  .listen(port);
