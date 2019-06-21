node_version = v10.16.0
node = node/bin/node
npm = node/bin/npm
tool_names = baiReadDepther bamstatsAlive craiReadDepther curl samtools
tools = $(patsubst %, tool_bin/%, $(tool_names))

export PATH := ./tools:$(PATH)

run: local_install
	$(node) index.js

local_install: node npm node_modules $(tools)

npm: node

node:
	wget https://nodejs.org/dist/$(node_version)/node-$(node_version)-linux-x64.tar.xz
	tar xf node-$(node_version)-linux-x64.tar.xz
	rm node-$(node_version)-linux-x64.tar.xz
	mv node-$(node_version)-linux-x64 node

tool_bin/%:
	./download_tool.sh $*

node_modules:
	$(npm) install

.PHONY: clean
clean:
	rm -rf node
	rm -rf node_modules
