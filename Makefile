node_version = v10.16.0
node = node/bin/node
npm = node/bin/npm
tool_names = baiReadDepther bamstatsAlive craiReadDepther samtools samtools_od bgzip tabix tabix_od vcfReadDepther coverage geneCoverage bcftools vt vep clinphen freebayes knownVariants_2
tools = $(patsubst %, tool_bin/%, $(tool_names))

export PATH := ./tools:$(PATH)

run: local_install
	$(node) src/index.js

deploy_aws: local_install
	./deploy_aws.sh

local_install: node npm node_modules tool_bin $(tools) static

npm: node

node:
	wget https://nodejs.org/dist/$(node_version)/node-$(node_version)-linux-x64.tar.xz
	tar xf node-$(node_version)-linux-x64.tar.xz
	rm node-$(node_version)-linux-x64.tar.xz
	mv node-$(node_version)-linux-x64 node

tool_bin:
	mkdir tool_bin

tool_bin/%:
	tools/download_tool.sh $*

static:
	./populate_static.sh

node_modules:
	$(node) $(npm) install

.PHONY: clean
clean:
	rm -rf node
	rm -rf node_modules
	rm -rf tool_bin
