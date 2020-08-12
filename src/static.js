const fs = require('fs');


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

  if (stats.isDirectory()) {
    ctx.status = 404;
    ctx.body = "Not Found";
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

  ctx.set('Accept-Ranges', 'bytes');
  ctx.set('Cache-Control', 'max-age=86400');

  ctx.body = stream;
}

module.exports = {
  serveStatic,
};
