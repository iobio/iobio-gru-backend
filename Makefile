node_version = v10.16.0
node = node/bin/node
npm = node/bin/npm

export PATH := ./tools:$(PATH)

run: local_install
	$(node) index.js

local_install: node npm node_modules

npm: node

node:
	wget https://nodejs.org/dist/$(node_version)/node-$(node_version)-linux-x64.tar.xz
	tar xf node-$(node_version)-linux-x64.tar.xz
	rm node-$(node_version)-linux-x64.tar.xz
	mv node-$(node_version)-linux-x64 node

node_modules:
	$(npm) install

.PHONY: clean
clean:
	rm -rf node
	rm -rf node_modules
