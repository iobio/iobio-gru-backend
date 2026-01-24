import fs from 'fs';
import path from 'path';
import mimelib from 'mime';
import { Readable } from 'node:stream';

async function serveStatic(req, fsPath) {

  const url = new URL(req.url);

  let stats;
  try {
    stats = await fs.promises.stat(fsPath);
  }
  catch (e) {
    return new Response("Not found", {
      status: 404,
    });
  }

  const lastModified = new Date(stats.mtimeMs).toUTCString();

  if (stats.isDirectory()) {

    let itemsHtml = '';
    const items = await fs.promises.readdir(fsPath, { withFileTypes: true });
    for (const item of items) {
      const name = item.isDirectory() ? item.name + '/' : item.name;
      itemsHtml += `  <div><a href='./${name}'>${name}</a></div>\n`;
    }

    const resHtml = `
      <!doctype html>
      <html>
      <head>
      </head>
      <body>
      ${itemsHtml}
      </body>
      </html>
    `

    return new Response(resHtml, {
      status: 200,
      headers: {
        'Content-Type': 'text/html',
      },
    });
  }

  const rangeHeader = req.headers.get('range');

  let status = 500;
  let headers = {};
  let nodeStream;

  // TODO: parse byte range specs properly according to
  // https://tools.ietf.org/html/rfc7233
  if (rangeHeader) {

    const range = {};
    const right = rangeHeader.split('=')[1];
    const rangeParts = right.split('-');
    range.start = Number(rangeParts[0]);
    range.end = stats.size - 1;

    if (rangeParts[1]) {
      range.end = Number(rangeParts[1]);
    }

    const originalSize = stats.size;

    nodeStream = fs.createReadStream(fsPath, {
      start: range.start,
      end: range.end,
    });

    status = 206;
    headers = {
      'Content-Range': `bytes ${range.start}-${range.end}/${originalSize}`,
      'Content-Length': range.end - range.start + 1,
    };
  }
  else {
    nodeStream = fs.createReadStream(fsPath);
  
    status = 200;
    headers = {
      'Content-Length': `${stats.size}`,
    };
  }

  const mime = mimelib.getType(path.extname(fsPath));
  if (mime) {
    headers['Content-Type'] = mime;
  }

  headers['Accept-Ranges'] = 'bytes';
  headers['Cache-Control'] = 'no-store';
  headers['Last-Modified'] = lastModified;

  let body = '';

  if (req.method === 'HEAD') {
    status = 200;
  }
  else {
    body = Readable.toWeb(nodeStream);
  }

  return new Response(body, {
    status,
    headers,
  });
}

export {
  serveStatic,
};
