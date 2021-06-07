node_version = v12.16.3

local_install: node node_modules

node:
	curl -LO https://nodejs.org/dist/$(node_version)/node-$(node_version)-linux-x64.tar.xz
	tar xf node-$(node_version)-linux-x64.tar.xz
	rm node-$(node_version)-linux-x64.tar.xz
	mv node-$(node_version)-linux-x64 node

# This forces the use of the local node. This is necessary because sqlite3 (and
# possibly other dependencies appears to be using the node in $PATH, even when
# calling `node/bin/node node/bin/npm install` directly. This leads to a
# version mismatch.
node_modules: export PATH := $(PWD)/node/bin:$(PATH)
node_modules: node
	npm install

.PHONY: clean
clean:
	rm -rf node
	rm -rf node_modules
