const fs = require('fs');
const path = require('path');
const { getType } = require('mime');


async function serveStatic(ctx, fsPath) {

  let stats;
  try {
    stats = await fs.promises.stat(fsPath);
  }
  catch (e) {
    ctx.status = 404;
    ctx.body = "Not Found";
    return;
  }

  const lastModified = new Date(stats.mtimeMs).toUTCString();

  if (stats.isDirectory()) {

    let itemsHtml = '';
    const items = await fs.promises.readdir(fsPath, { withFileTypes: true });
    for (const item of items) {
      const name = item.isDirectory() ? item.name + '/' : item.name;
      itemsHtml += `  <div><a href='./${name}'>${name}</a></div>\n`;
    }

    resHtml = `
      <!doctype html>
      <html>
      <head>
      </head>
      <body>
      ${itemsHtml}
      </body>
      </html>
    `

    ctx.set('Content-Type', 'text/html');

    ctx.status = 200;
    ctx.body = resHtml;
    return;
  }

  const rangeHeader = ctx.headers['range'];

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

    ctx.set('Content-Range', `bytes ${range.start}-${range.end}/${originalSize}`);
    ctx.set('Content-Length', range.end - range.start + 1);
    ctx.status = 206;

    stream = fs.createReadStream(fsPath, {
      start: range.start,
      end: range.end,
    });
  }
  else {
    ctx.set('Content-Length', `${stats.size}`);
    stream = fs.createReadStream(fsPath);
  }

  const mime = getType(path.extname(fsPath));
  if (mime) {
    ctx.set('Content-Type', mime);
  }

  ctx.set('Accept-Ranges', 'bytes');
  ctx.set('Cache-Control', 'no-store');
  ctx.set('Last-Modified', lastModified);

  if (ctx.method == 'HEAD') {
    ctx.status = 200;
  }
  else {
    ctx.body = stream;
  }
}

module.exports = {
  serveStatic,
};
