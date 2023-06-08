const os = require('os');
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const bodyParser = require('koa-bodyparser');
const path = require('path');
const spawn = require('child_process').spawn;
const process = require('process');
const gene2PhenoRouter = require('./gene2pheno.js');
const geneInfoRouter = require('./geneinfo.js');
const genomeBuildRouter = require('./genomebuild.js');
const hpoRouter = require('./hpo.js');
const { parseArgs, dataPath } = require('./utils.js');
const fs = require('fs');
const { serveStatic } = require('./static.js');
const stream = require('stream');
const semver = require('semver');

const MAX_STDERR_LEN = 1048576;
const MIN_DATA_DIR_VERSION = '1.9.0';

console.log(`Using data directory ${path.resolve(dataPath(''))}`);
const dataDirVersion = fs.readFileSync(dataPath('VERSION')).toString();
if (semver.lt(dataDirVersion, MIN_DATA_DIR_VERSION)) {
  console.error(`Data directory must be at least version ${MIN_DATA_DIR_VERSION} (found ${dataDirVersion})`);
  process.exit(1);
}


// Clean up any older tmp files that may have been left behind after a crash
const tmpFiles = fs.readdirSync(os.tmpdir());
for (const entName of tmpFiles) {
  if (entName.startsWith('gru-')) {
    const tmpPath = path.join(os.tmpdir(), entName);
    console.log(`Cleaning up tmp file ${tmpPath}`);
    fs.rmSync(tmpPath, { recursive: true, force: true });
  }
}

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
    const indexUrl = params.indexUrl ? params.indexUrl : '';
    await handle(ctx, 'getChromosomes.sh', [params.url, indexUrl]);
});

router.post('/vcfReadDepth', async (ctx) => {
    const params = JSON.parse(ctx.request.body);
    await handle(ctx, 'vcfReadDepth.sh', [params.url]);
});

router.post('/alignmentCoverage', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

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

    await handle(ctx, 'getClinvarVariants.sh', args, { ignoreStderr: true });
});

router.post('/getClinvarVariantsV2', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
    const contigStr = genContigFileStr(params.refNames);
    const regionStr = genRegionsStr(params.regions);
    const refFastaFile = dataPath(params.refFastaFile);
    const gnomadMergeAnnots = params.gnomadMergeAnnots ? params.gnomadMergeAnnots : '';

    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr, refFastaFile, 
        params.genomeBuildName, gnomadMergeAnnots, params.clinSigFilterPhrase
    ];

    await handle(ctx, 'getClinvarVariantsV2.sh', args, { ignoreStderr: true });
});

router.post('/annotateVariants', async (ctx) => {

    const params = JSON.parse(ctx.request.body);

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

    await handle(ctx, 'annotateVariants.sh', args, { ignoreStderr: true });
});


router.post('/annotateVariantsV2', async (ctx) => {

    const params = JSON.parse(ctx.request.body);

    const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
    const contigStr = genContigFileStr(params.refNames);
    const regionStr = genRegionsStr(params.regions);
    const vcfSampleNamesStr = params.vcfSampleNames.join("\n");
    const refFastaFile = dataPath(params.refFastaFile);
    const vepCacheDir = dataPath('vep-cache');
    const vepREVELFile = dataPath(params.vepREVELFile);
    const vepPluginDir = dataPath('vep-cache/Plugins');
    const gnomadMergeAnnots = params.gnomadMergeAnnots ? params.gnomadMergeAnnots : '';

    const args = [
        params.vcfUrl, tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
        refFastaFile, params.genomeBuildName, vepCacheDir, vepREVELFile, params.vepAF,
        vepPluginDir, params.hgvsNotation, params.getRsId, gnomadMergeAnnots, 
        params.decompose

    ];

    await handle(ctx, 'annotateVariantsV2.sh', args, { ignoreStderr: true });
});
router.post('/annotateEnrichmentCounts', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

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

// secondary option used by oncogene
router.post('/annotateSomaticVariantsVep', async (ctx) => {
  const params = JSON.parse(ctx.request.body);
  const vepCacheDir = dataPath('vep-cache');

  let refFastaFile = dataPath('references/GRCh37/human_g1k_v37_decoy_phix.fasta');
  if (params.genomeBuildName === 'GRCh38') {
    refFastaFile = dataPath('references/GRCh38/human_g1k_v38_decoy_phix.fasta');
  }
  
  const args = [params.vcfUrl, params.selectedSamplesStr, params.geneRegionsStr, params.somaticFilterPhrase, params.genomeBuildName, vepCacheDir, refFastaFile];
  
  await handle(ctx, 'annotateSomaticVariantsVep.sh', args, { ignoreStderr: true });
});

// primary, faster option used by oncogene
router.post('/annotateSomaticVariantsBcsq', async (ctx) => {
  const params = JSON.parse(ctx.request.body);

  let refFastaFile = dataPath('references/GRCh37/human_g1k_v37_decoy_phix.fasta');
  let gffFile = '/iobio-gru-backend/static/ensembl/GRCh37/geneSet37.gff3.gz';

  if (params.genomeBuildName === 'GRCh38') {
    refFastaFile = dataPath('references/GRCh38/human_g1k_v38_decoy_phix.fasta');
    gffFile = '/iobio-gru-backend/static/ensembl/GRCh38/geneSet38.gff3.gz';
  }

  const args = [params.vcfUrl, params.selectedSamplesStr, params.geneRegionsStr, params.somaticFilterPhrase, refFastaFile, gffFile];
  
  await handle(ctx, 'annotateSomaticVariantsBcsq.sh', args, { ignoreStderr: true });
});

router.post('/freebayesJointCall', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

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

  const extraArgs = freebayesArgs;

  const args = [
    alignments, indices, region, refFastaFile, useSuggestedVariants,
    params.clinvarUrl, params.genomeBuildName, vepREVELFile, params.vepAF,
    samplesFileStr, extraArgs, vepCacheDir, vepPluginDir,
    gnomadUrl, gnomadRegionStr, gnomadHeaderFile, decompose, dataPath(''),
  ];

  await handle(ctx, 'freebayesJointCall.sh', args, { ignoreStderr: true });
});


router.post('/freebayesJointCallV2', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

  const alignments = params.alignmentSources.map(aln => aln.bamUrl).join(',');
  const indices = params.alignmentSources.map(aln => aln.baiUrl).join(',');
  const region = genRegionStr(params.region);
  const vepREVELFile = dataPath(params.vepREVELFile);
  const refFastaFile = dataPath(params.refFastaFile);
  const contigStr = genContigFileStr(params.refNames);
  const samplesFileStr = params.sampleNames.join('\n');

  const vepCacheDir = dataPath('vep-cache');
  const vepPluginDir = dataPath('vep-cache/Plugins');

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

  const extraArgs = freebayesArgs;

  const dataDir = dataPath('');

  const args = [
    alignments, indices, region, refFastaFile, useSuggestedVariants,
    params.clinvarUrl, params.genomeBuildName, vepREVELFile, params.vepAF,
    samplesFileStr, extraArgs, vepCacheDir, vepPluginDir,
    decompose, contigStr, dataDir
  ];

  await handle(ctx, 'freebayesJointCallV2.sh', args, { ignoreStderr: true });
});


router.post('/clinvarCountsForGene', async (ctx) => {
  const params = JSON.parse(ctx.request.body);

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
  const tmpDir = await fs.promises.mkdtemp(path.join(os.tmpdir(), 'gru-'));
  const tmpFilePath = path.join(tmpDir, 'clin_report');
  await fs.promises.writeFile(tmpFilePath, ctx.request.body);
  const args = [tmpFilePath];
  await handle(ctx, 'clinReport.sh', args); 	
  await fs.promises.rm(tmpDir, { recursive: true, force: true });
});
	



// oncogene endpoints
//
router.post('/getIdColumns', async (ctx) => {

    const params = JSON.parse(ctx.request.body);

    const regionStr = genRegionsStr(params.regions, ",");
    const args = [
        params.vcfUrl, regionStr
    ];

    await handle(ctx, 'getIdColumns.sh', args, { ignoreStderr: true });
});

router.post('/checkBamBai', async (ctx) => {
    const params = JSON.parse(ctx.request.body);

    const indexUrl = params.indexUrl ? params.indexUrl : '';

    const args = [ params.url, indexUrl, params.region, dataPath('') ];
    await handle(ctx, 'checkBamBai.sh', args, { ignoreStderr: true });
});



// vcf.iobio endpoints
router.post('/vcfStatsStream', async (ctx) => {

  const params = JSON.parse(ctx.request.body);

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

  await handle(ctx, 'vcfStatsStream.sh', args, { ignoreStderr: true });
}); 


async function handle(ctx, scriptName, args, options) {

  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'gru-'));

  if (ctx.gruParams) {
    fs.writeFileSync(path.join(tmpDir, `gru_params_${scriptName}.json`), JSON.stringify(ctx.gruParams, null, 2));
  }

  const opts = {cwd: tmpDir, ...options};
  
  const scriptPath = path.join(__dirname, '../scripts', scriptName);
  const proc = spawn(scriptPath, args, opts);

  // Kill process if it runs for more than 5 minutes
  const MINUTE_MS = 60*1000;
  const timeoutId = setTimeout(() => {
    console.error("Timed out. Killing process for request", ctx.gruParams._requestId);
    proc.kill('SIGKILL');
  }, 5 * MINUTE_MS);

  const out = stream.PassThrough();

  ctx.body = out;

  let closed = false;

  ctx.res.once('close', () => {
    closed = true;
    proc.stdout.destroy();
  });

  proc.stdout.on('data', (chunk) => {
    if (!closed) {
      out.write(chunk);
    }
  });

  let stderr = "";
  proc.stderr.on('data', (chunk) => {
    if (stderr.length < MAX_STDERR_LEN) {
      stderr += chunk;
    }
  });

  return new Promise((resolve, reject) => {
    proc.on('exit', (exitCode) => {

      clearTimeout(timeoutId);

      fs.rmSync(tmpDir, { recursive: true, force: true });

      if (exitCode !== 0) {
        const timestamp = new Date().toISOString();
        console.log(`${timestamp}\t${ctx.gruParams._requestId}\terror\t${ctx.url}`);
        console.log("stderr:");
        console.log(stderr);
        console.log("params:");
        console.log(ctx.gruParams);

        if (ctx.gruParams._appendErrors === true) {
          out.write("GRU_ERROR_SENTINEL");
          out.write(JSON.stringify({
            stderr,
          }));
        }
      }

      out.end();
      resolve();
    });
  });
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

async function logger(ctx, next) {
  const contentType = ctx.get('content-type').split(';')[0];

  if (ctx.method !== 'POST' || contentType != 'text/plain') {
    await next();
    return;
  }

  const params = JSON.parse(ctx.request.body);

  ctx.gruParams = params;

  let timestamp = new Date().toISOString();
  const start = Date.now();
  console.log(`${timestamp}\t${params._requestId}\tstart\t${ctx.url}\t${params._attemptNum}`);

  await next();

  // Modifed from https://github.com/koajs/logger/blob/f8edfa00cb5af7e696cf276ddc8b482accd9f7a9/index.js#L76
  const { res } = ctx;

  const onFinish = done.bind(null, 'finish');
  const onClose = done.bind(null, 'close');

  res.once('finish', onFinish);
  res.once('close', onClose);

  function done(evt) {
    res.removeListener('finish', onFinish);
    res.removeListener('close', onClose);

    const message = evt === 'finish' ? 'finish' : 'canceled';

    timestamp = new Date().toISOString();
    const seconds = (Date.now() - start) / 1000;
    console.log(`${timestamp}\t${params._requestId}\t${message}\t${ctx.url}\t${seconds} seconds`);
  }
}

const server = new Koa();
server
  .use(cors({
    origin: '*',
    maxAge: 86400,
  }))
  .use(bodyParser({
    enableTypes: ['json', 'text'],
    jsonLimit: '10mb',
    textLimit: '10mb',
  }))
  .use(logger)
  .use(rootRouter.routes())
  .use(rootRouter.allowedMethods())
  .listen(port);
