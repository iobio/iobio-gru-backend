import { spawn } from 'node:child_process';
import fs from 'fs/promises';
import path from 'node:path';
import os from 'node:os';

async function phenolyzerHandler(req, url = new URL(req.url)) {
  const params = {
    term: url.searchParams.get("term"),
    refresh: url.searchParams.get("refresh"),
  };

  if (!params.term) {
    return new Response("Must provide term parameter", {
      status: 400,
    });
  }

  const cachePath = buildCachePath(params.term);

  const { readable, writable } = new TransformStream();
  const writer = writable.getWriter();

  const headers = {
    'Cache-Control': 'max-age=86400',
    'Content-Type': 'application/json',
  };

  if (params.refresh === 'true') {
    try {
      await fs.unlink(cachePath);
    }
    catch (e) {
    }
  }

  let stats;
  try {
    stats = await fs.stat(cachePath);
  }
  catch(e) {
  }

  if (stats && stats.size === 0) {
    if (pendingTooLong(stats)) {
      try {
        await fs.unlink(cachePath);
      }
      catch (e) {
      }

      stats = null;
    }
    else {
      const body = JSON.stringify({
        record: 'pending',
      });

      return new Response(body, {
        headers: {
          'Content-Type': 'application/json',
          'Cache-Control': 'no-store',
        },
      });
    }
  } 

  if (stats) {
    const data = await fs.readFile(cachePath, 'utf8');

    const body = JSON.stringify({
      record: data,
    });

    return new Response(body, {
      headers,
    });
  }
  else {

    await ensureCacheDir(params.term);
    await fs.writeFile(cachePath, '');

    const tmpDir = await fs.mkdtemp(path.join(os.tmpdir(), 'gru-'));

    const proc = spawn('phenolyzer', [params.term], { cwd: tmpDir });
    proc.stdout.setEncoding('utf8');

    let data = '';
    proc.stdout.on('data', (chunk) => {
      data += chunk;
    });

    let ended = false;
    proc.stdout.on('end', async () => {
      ended = true;
    });

    proc.on('exit', async () => {

      fs.rm(tmpDir, { recursive: true, force: true });

      if (!ended) {
        console.error("Attempted to write before stream ended");
      }

      if (proc.exitCode === 0) {
        await fs.writeFile(cachePath, data);
      }
      else {
        console.error("Phenolyzer failed for term:", `"${params.term}"`);
      }
    });

    const body = JSON.stringify({
      record: 'queued',
    });

    return new Response(body, {
      headers: {
        'Content-Type': 'application/json',
        'Cache-Control': 'no-store',
      },
    });
  }

  return new Response("Server error", { status: 500 });
}

function pendingTooLong(stats) {
  const fiveMinutes = 5 * 60 * 1000;
  const now = new Date();
  return (now - stats.mtime.getTime()) > fiveMinutes;
}

async function ensureCacheDir(term) {
  const encTerm = encodeURIComponent(term);
  const cacheDirname = encTerm.slice(0, 2);
  const cacheDir = path.join('cache', cacheDirname);
  try {
    await fs.mkdir(cacheDir, { recursive: true });
  }
  catch (e) {
    // no-op
  }
}

function buildCachePath(term) {
  const encTerm = encodeURIComponent(term);
  const cacheDirname = encTerm.slice(0, 2);
  const cacheDir = path.join('cache', cacheDirname);
  const cachePath = path.join(cacheDir, encTerm);
  return cachePath;
}

export {
  phenolyzerHandler,
};
