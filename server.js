const { RPCServer, call } = require('./src/index');

const server = new RPCServer({});
server.start(9001);
