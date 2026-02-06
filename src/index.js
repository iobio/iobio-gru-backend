import os from 'node:os';
import path from 'node:path';
import { spawn } from 'node:child_process';
import process from 'node:process';
import { gene2PhenoHandler } from './gene2pheno.js';
import { createGeneInfoHandler } from './geneinfo.js';
import { genomeBuildHandler } from './genomebuild.js';
import { hpoHandler } from './hpo.js';
import { preciseReadDepthHandler } from './precise_read_depth_handler.js';
import { parseArgs, dataPath, replaceAll, genRegionStr, semverLess } from './utils.js';
import fs from 'node:fs';
import { serveStatic } from './static.js';
import { fileURLToPath } from 'node:url';
//import AutoEncrypt from '@small-tech/auto-encrypt';
import { serve } from '@anderspitman/fetch-handler';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const MAX_STDERR_LEN = 1048576;
const MIN_DATA_DIR_VERSION = '1.12.0';

const gruVersion = fs.readFileSync(path.join(__dirname, '..', 'VERSION')).toString().trim();

console.log(`Using data directory ${path.resolve(dataPath(''))}`);
const dataDirVersion = fs.readFileSync(dataPath('VERSION')).toString().trim();
if (semverLess(dataDirVersion, MIN_DATA_DIR_VERSION)) {
  console.error(`Data directory must be at least version ${MIN_DATA_DIR_VERSION} (found ${dataDirVersion})`);
  process.exit(1);
}

const statusData = {
  service_description: "iobio gru backend server",
  version: gruVersion,
  data_version: dataDirVersion,
};


// Clean up any older tmp files that may have been left behind after a crash
const tmpFiles = fs.readdirSync(os.tmpdir());
for (const entName of tmpFiles) {
  if (entName.startsWith('gru-')) {
    const tmpPath = path.join(os.tmpdir(), entName);
    console.log(`Cleaning up tmp file ${tmpPath}`);
    fs.rmSync(tmpPath, { recursive: true, force: true });
  }
}

// TODO: port and test embedded app functionality
//router.get('/static/*', async (ctx) => {
//
//  // TODO: Hack to remove prefix when hosting app within gru because I don't
//  // really understand how to do this properly with koa-router
//  const reqPath = ctx.path.startsWith('/gru') ? ctx.path.slice(4) : ctx.path;
//
//  const fsPath = path.join(__dirname, '..', reqPath);
//  await serveStatic(ctx, fsPath);
//  ctx.set('Cache-Control', 'max-age=86400');
//});


function runScript(req, scriptName, args, gruParams, options) {

  const url = new URL(req.url);

  const scriptPath = path.join(__dirname, '../scripts', scriptName);

  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'gru-'));

  const opts = {cwd: tmpDir, ...options};
  
  const proc = spawn(scriptPath, args, opts);

  // Kill process if it runs for more than 5 minutes
  const MINUTE_MS = 60*1000;
  const timeoutId = setTimeout(() => {
    console.error("Timed out. Killing process for request", gruParams._requestId);
    proc.kill('SIGKILL');
  }, 5 * MINUTE_MS);

  const { readable, writable } = new TransformStream();

  const writer = writable.getWriter();

  let closed = false;

  req.signal.addEventListener('abort', (evt) => {
    closed = true;
    proc.stdout.destroy();
  });

  proc.stdout.on('data', async (chunk) => {
    if (!closed) {
      await writer.write(chunk);
    }
  });

  let stderr = "";
  proc.stderr.on('data', (chunk) => {
    if (stderr.length < MAX_STDERR_LEN) {
      stderr += chunk;
    }
  });

  proc.on('exit', async (exitCode) => {

    clearTimeout(timeoutId);

    fs.rmSync(tmpDir, { recursive: true, force: true });

    if (exitCode !== 0) {
      const timestamp = new Date().toISOString();
      console.log(`${timestamp}\t${gruParams._requestId}\terror\t${req.url}`);
      console.log("stderr:");
      console.log(stderr);
      console.log("params:");
      console.log(gruParams);

      if (gruParams._appendErrors === true) {
        if (!writer.closed) {
          await writer.write("GRU_ERROR_SENTINEL");
          await writer.write(JSON.stringify({
            stderr,
          }));
        }
      }
    }

    if (!closed) {
      writer.close();
    }
  });

  const res = new Response(readable);
  res.headers.set('Access-Control-Allow-Origin', '*');

  return res;
}

function genContigFileStr(refNames) {
  let contigStr = "";
  for (const ref of refNames) {
    contigStr += "##contig=<ID=" + ref + ">\n";
  }
  return contigStr;
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

const enableHttps = args['--enable-https'] === 'true';


// TODO: port and test embedded app functionality
//import Router from 'koa-router';
//const router = new Router();
//let rootRouter = router;
//
// Allows a frontend app to be hosted in the same process. This is particularly
// useful if you want to run an entire app frontend and backend in a single
// docker container.
//if (args['--app-dir']) {
//  rootRouter = new Router();
//  rootRouter.use('/gru', router.routes(), router.allowedMethods());
//
//  rootRouter.get('/gru', async (ctx) => {
//    ctx.body = statusData;
//  });
//
//  rootRouter.get('/*', async (ctx, next) => {
//
//    const fsPath = path.join(args['--app-dir'], ctx.path);
//    await serveStatic(ctx, fsPath);
//
//    // Bit of a hack. If previous attempt to serve didn't find the file,
//    // default to the root index.html.
//    if (ctx.status === 404) {
//      const fsPath = path.join(args['--app-dir'], 'index.html');
//      await serveStatic(ctx, fsPath);
//    }
//  });
//}
//else {
//  rootRouter.get('/', async (ctx) => {
//    ctx.body = statusData;
//  });
//}

function wrapper(handler) {
  return async (req, nodeReq, nodeRes) => {

    const url = new URL(req.url);

    const contentType = req.headers.get('content-type')?.split(';')[0];

    let timestamp = new Date().toISOString();

    if (req.method !== 'POST' || contentType != 'text/plain') {
      console.log(`${timestamp}\t${req.method}\t${url.pathname}`);
      const res = await handler(req, nodeReq, nodeRes);
      return res;
    }

    const params = new URLSearchParams(url.search);

    const start = Date.now();
    console.log(`${timestamp}\t${params._requestId}\tstart\t${url.pathname}\t${params._attemptNum}`);

    const res = await handler(req, nodeReq, nodeRes);

    if (res && res.body) {
      return withCompletedCallback(res, () => {
        timestamp = new Date().toISOString();
        const seconds = (Date.now() - start) / 1000;
        console.log(`${timestamp}\t${params._requestId}\tfinish\t${url.pathname}\t${seconds} seconds`);
      });
    }
    else {
      return res;
    }
  };
}


// TODO: figure out HTTPS
//if (enableHttps) {
//
//  const domain = args['--domain'];
//
//  const serverOptions = {
//    domains: [domain],
//  };
//
//  const server = AutoEncrypt.https.createServer(serverOptions, app.callback());
//
//  server.listen(koaPort);
//}
//else {
//  app.listen(koaPort);
//}


function addCors(handler) {
  return async (req) => {
    const res = await handler(req);
    res.headers.set('Access-Control-Allow-Origin', '*');
    return res;
  };
}

const geneInfoHandler = addCors(createGeneInfoHandler({ prefix: '/geneinfo' }));
const g2pHandler = addCors(gene2PhenoHandler);

function withCompletedCallback(res, callback) {
  if (!res.body) {
    callback();
    return res;
  }

  const { readable, writable } = new TransformStream();

  res.body.pipeTo(writable)
    .catch((e) => {
      console.error("Error in withCompletedCallback stream:", e);
    })
    .finally(() => {
      callback();
    });

  return new Response(readable, res);
}

async function handler(req, nodeReq, nodeRes) {

  const url = new URL(req.url);

  if (url.pathname === '/') {
    return new Response(JSON.stringify(statusData), {
      headers: {
        'Content-Type': 'application/json',
      },
    });
  }
  else if (url.pathname.startsWith('/static')) {
    const fsPath = path.join(__dirname, '..', url.pathname);
    return serveStatic(req, fsPath);
  }
  else if (url.pathname.startsWith('/gene2pheno')) {
    return g2pHandler(req);
  }
  else if (url.pathname.startsWith('/geneinfo')) {
    return geneInfoHandler(req);
  }
  else if (url.pathname.startsWith('/genomebuild')) {
    return genomeBuildHandler(req);
  }
  else if (url.pathname.startsWith('/hpo/hot/lookup')) {
    return hpoHandler(req);
  }
  else if (url.pathname.startsWith('/gru_data/')) {
    const fsPath = path.join(dataPath(''), url.pathname.slice(10));
    return serveStatic(req, fsPath);
  }
  else if (url.pathname === '/clinphen') {
    const notes = url.searchParams.get('notes');
    return runScript(req, 'clinphen.sh', [notes], { notes });
  }
  else if (url.pathname === '/phenotypeExtractor') {
    const notes = url.searchParams.get('notes');
    return runScript(req, 'phenotypeExtractor.sh', [notes], { notes });
  }
  else if (url.pathname === '/clinReport') {
    // Copy the data into a temporary file and then pass the path. It was failing
    // before, I'm pretty sure because the file was too large (~3MB) to pass
    // through the child spawing interface.
    const body = await req.text();
    const gruParams = { _requestId: new URL(req.url).searchParams.get('_requestId') };
    const tmpDir = await fs.promises.mkdtemp(path.join(os.tmpdir(), 'gru-'));
    const tmpFilePath = path.join(tmpDir, 'clin_report');
    await fs.promises.writeFile(tmpFilePath, body);
    const args = [tmpFilePath];
    const res = runScript(req, 'clinReport.sh', args, gruParams);
    return withCompletedCallback(res, () => {
      fs.promises.rm(tmpDir, { recursive: true, force: true });
    });
  }
  else {
    let params;

    try {
      params = await req.json();
    }
    catch (e) {
      return new Response("Invalid JSON payload", {
        status: 400,
      });
    }

    if (url.pathname === '/viewAlignments') {
      const args = [params.url];

      if (params.regions) {
        const samtoolsRegions = genRegionsStr(params.regions);
        args.push(samtoolsRegions);
      }
      return runScript(req, 'viewAlignments.sh', args, params);
    }
    else if (url.pathname === '/alignmentHeader') {
      return runScript(req, 'alignmentHeader.sh', [params.url], params);
    }
    else if (url.pathname === '/baiReadDepth') {
        return runScript(req, 'baiReadDepth.sh', [params.url], params);
    }
    else if (url.pathname === '/craiReadDepth') {
        return runScript(req, 'craiReadDepth.sh', [params.url], params);
    }
    else if (url.pathname === '/alignmentStatsStream') {
        const samtoolsRegions = genRegionsStr(params.regions);
        const bamstatsRegions = JSON.stringify(params.regions.map(function(d) { return {start:d.start,end:d.end,chr:d.name};}));
        const indexUrl = params.indexUrl ? params.indexUrl : '';
        const args = [
            params.url,
            indexUrl,
            samtoolsRegions,
            bamstatsRegions,
            dataPath(''),
        ];
        return runScript(req, 'alignmentStatsStream.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/bedRegion') {
        const url = params.url;
        const indexUrl = params.indexUrl ? params.indexUrl : '';
        const region = genRegionStr(params.region);
        return runScript(req, 'bedRegion.sh', [params.url, indexUrl, region], params);
    }
    else if (url.pathname === '/bigWigDepther') {
        const url = params.url;
        const region = genRegionStr(params.region);
        return runScript(req, 'bigWigDepther.sh', ["--url", params.url, "--region", region], params);
    }
    else if (url.pathname === '/getReferenceSequence') {
        const refFastaPath = dataPath(params.fastaPath);
        const region = genRegionStr(params.region);
        return runScript(req, 'getReferenceSequence.sh', [refFastaPath, region], params);
    }
    else if (url.pathname === '/variantHeader') {
      const indexUrl = params.indexUrl ? params.indexUrl : '';
      return runScript(req, 'variantHeader.sh', [params.url, indexUrl], params);
    }
    else if (url.pathname === '/getChromosomes') {
        const indexUrl = params.indexUrl ? params.indexUrl : '';
        return runScript(req, 'getChromosomes.sh', [params.url, indexUrl], params);
    }
    else if (url.pathname === '/vcfReadDepth') {
        return runScript(req, 'vcfReadDepth.sh', [params.url], params);
    }
    else if (url.pathname === '/alignmentCoverage') {
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
        return runScript(req, 'alignmentCoverage.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/geneCoverage') {
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
        return runScript(req, 'geneCoverage.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/normalizeVariants') {
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
        return runScript(req, 'normalizeVariants.sh', args, params);
    }
    else if (url.pathname === '/getClinvarVariants') {
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

        return runScript(req, 'getClinvarVariants.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/getClinvarVariantsV2') {
        const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
        const contigStr = genContigFileStr(params.refNames);
        const regionStr = genRegionsStr(params.regions);
        const refFastaFile = dataPath(params.refFastaFile);
        const gnomadMergeAnnots = params.gnomadMergeAnnots ? params.gnomadMergeAnnots : '';

        const args = [
            params.vcfUrl, tbiUrl, regionStr, contigStr, refFastaFile,
            params.genomeBuildName, gnomadMergeAnnots, params.clinSigFilterPhrase
        ];

        return runScript(req, 'getClinvarVariantsV2.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateVariants') {
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
        return runScript(req, 'annotateVariants.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateVariantsV2') {
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
        return runScript(req, 'annotateVariantsV2.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateVariantsV3') {
        const tbiUrl = params.tbiUrl ? params.tbiUrl : '';
        const contigStr = genContigFileStr(params.refNames);
        const regionStr = genRegionsStr(params.regions);
        const vcfSampleNamesStr = params.vcfSampleNames.join("\n");
        const refFastaFile = dataPath(params.refFastaFile);
        const vepCacheDir = dataPath('vep-cache');
        const vepREVELFile = dataPath(params.vepREVELFile);
        const vepPluginDir = dataPath('vep-cache/Plugins');

        const args = [
            params.vcfUrl, tbiUrl, regionStr, contigStr, vcfSampleNamesStr,
            refFastaFile, params.genomeBuildName, vepCacheDir, vepREVELFile,
            vepPluginDir, params.hgvsNotation, params.getRsId, params.decompose
        ];
        return runScript(req, 'annotateVariantsV3.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateEnrichmentCounts') {
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
        return runScript(req, 'annotateEnrichmentCounts.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateSomaticVariantsVep') {
        const vepCacheDir = dataPath('vep-cache');
        let refFastaFile = dataPath('references/GRCh37/human_g1k_v37_decoy_phix.fasta');
        if (params.genomeBuildName === 'GRCh38') {
            refFastaFile = dataPath('references/GRCh38/human_g1k_v38_decoy_phix.fasta');
        }
        const args = [params.vcfUrl, params.selectedSamplesStr, params.geneRegionsStr, params.somaticFilterPhrase, params.genomeBuildName, vepCacheDir, refFastaFile];
        return runScript(req, 'annotateSomaticVariantsVep.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/annotateSomaticVariantsBcsq') {
        let refFastaFile = dataPath('references/GRCh37/human_g1k_v37_decoy_phix.fasta');
        let gffFile = '/iobio-gru-backend/static/ensembl/GRCh37/geneSet37.gff3.gz';

        if (params.genomeBuildName === 'GRCh38') {
            refFastaFile = dataPath('references/GRCh38/human_g1k_v38_decoy_phix.fasta');
            gffFile = '/iobio-gru-backend/static/ensembl/GRCh38/geneSet38.gff3.gz';
        }
        const args = [params.vcfUrl, params.selectedSamplesStr, params.geneRegionsStr, params.somaticFilterPhrase, refFastaFile, gffFile];
        return runScript(req, 'annotateSomaticVariantsBcsq.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/freebayesJointCall') {
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
        return runScript(req, 'freebayesJointCall.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/freebayesJointCallV2') {
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
        return runScript(req, 'freebayesJointCallV2.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/freebayesJointCallV3') {
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
            params.clinvarUrl, params.genomeBuildName, vepREVELFile,
            samplesFileStr, extraArgs, vepCacheDir, vepPluginDir,
            decompose, contigStr, dataDir
        ];
        return runScript(req, 'freebayesJointCallV3.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/clinvarCountsForGene') {
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

        return runScript(req, 'clinvarCountsForGene.sh', args, params);
    }
    else if (url.pathname === '/bcftoolsView') {
        const regionStr = params.regions == undefined ? "" : genRegionsStr(params.regions, ",");
        const numRecords = params.numRecords == undefined ? "" : params.numRecords;
        const args = [
            params.vcfUrl, regionStr, numRecords
        ];
        return runScript(req, 'bcftoolsView.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/checkBamBai') {
        const indexUrl = params.indexUrl ? params.indexUrl : '';
        const args = [ params.url, indexUrl, params.region, dataPath('') ];
        return runScript(req, 'checkBamBai.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/vcfStatsStream') {
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
        return runScript(req, 'vcfStatsStream.sh', args, params, { ignoreStderr: true });
    }
    else if (url.pathname === '/preciseReadDepth') {
      return preciseReadDepthHandler(req, params);
    }
  }

  return new Response("Invalid request", {
    status: 400,
  });
}

serve({
  handler: wrapper(handler),
  //handler,
  port,
});
