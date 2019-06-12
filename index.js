const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const path = require('path');
const { run } = require('./process.js');

const spawn = require('child_process').spawn;
const process = require('process');


const app = new Koa();
const router = new Router();

router.get('/indexReadDepth', async (ctx) => {
  try {
    const scriptPath = path.join(__dirname, '/scripts/indexReadDepth.sh');
    const proc = await run(scriptPath, [ctx.query.url]);
    ctx.body = proc.stdout;
  }
  catch (e) {
    console.error(e);
    ctx.status = 500;
  }
});


app
  .use(cors())
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(9001);
