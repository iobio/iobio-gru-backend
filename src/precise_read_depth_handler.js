import path from 'node:path';
import { spawn } from 'node:child_process';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';

import { genRegionStr, dataPath } from './utils.js';


const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);


function preciseReadDepthHandler(req, params) {

  // +1 because samtools regions are 1-based
  const numBases = 1 + params.region.end - params.region.start;

  if (numBases > 20000000) {
    throw new Error("region too large");
  }

  if (numBases < params.numBins) {
    throw new Error("end-start must be more than numBins");
  }

  const samtoolsRegion = genRegionStr(params.region);

  const indexUrl = params.indexUrl ? params.indexUrl : '';

  const scriptPath = path.join(__dirname, '../scripts', 'preciseReadDepth.sh');
  const args = [
    params.url,
    indexUrl,
    samtoolsRegion,
    dataPath(''),
  ];
  const proc = spawn(scriptPath, args, {});

  const binSize = Math.ceil(numBases / params.numBins);

  const { readable, writable } = new TransformStream();

  const writer = writable.getWriter();

  const rl = readline.createInterface({
    input: proc.stdout,
  });

  let aggDepth = 0;
  let binIdx = 0;
  rl.on('line', async (line) => {
    const parts = line.split('\t');
    const base = parseInt(parts[1]);
    const depth = parseInt(parts[2]);
    aggDepth += depth;

    const baseIdx = base - params.region.start;
    const assignedBin = Math.floor((baseIdx / numBases) * params.numBins);

    if (assignedBin > binIdx || binIdx === params.numBins - 1) {
      binIdx += 1;

      const avgDepth = Math.floor(aggDepth / binSize);
      await writer.write(`${avgDepth}\n`);

      aggDepth = 0;
    }
  });

  rl.on('close', async () => {
    await writer.close();
  });

  return new Response(readable);
}

export {
  preciseReadDepthHandler,
};
