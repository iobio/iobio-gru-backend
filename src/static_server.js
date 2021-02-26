#!/usr/bin/env node

const path = require('path');
const Koa = require('koa');
const Router = require('koa-router');
const cors = require('@koa/cors');
const { serveStatic } = require('./static.js');

const router = new Router();

router.get('/*', async (ctx) => {
  const fsPath = path.join('.', ctx.path);
  await serveStatic(ctx, fsPath);
});

const port = process.argv[2] ? process.argv[2] : 8080;

const server = new Koa();
server
  .use(cors({
    origin: '*',
    maxAge: 86400,
  }))
  .use(router.routes())
  .use(router.allowedMethods())
  .listen(port);
