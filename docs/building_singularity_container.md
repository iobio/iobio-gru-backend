This document will outline the steps to build a singularity container. 

**Building a singularity container requires the following:** 
* A Definition File 
* Singularity build command 




**Creating a Definition file:** 
A singularity definition file (.def) is a set of blueprints explaining how to build a custom container. Follow the procedure to create a .def file.
1. Go to tools directory. example: `cd adit/iobio-gru-backend/tools` 
2. Create a directory for your tool. Let us assume that the name of our tool is clin service: `mkdir clinService`
3. Go to the newly created director: `cd clinService`
4. Create a def file: `vi clinService.def`
5. You can use this as a template to create a .def file: 
```
Bootstrap: library

From: alpine:3.9

%post

	apk update

	apk add git nodejs npm 

	apk add python python-dev build-base 

	apk add py-pip

	pip install nltk

	python -m nltk.downloader -d /usr/share/nltk_data wordnet

	git clone https://github.com/adityaekawade/clinService.git

	git checkout *commitId*

	cd clinService

	npm install

%runscript

	/Phenotype-extractor/cli.js $@
```
The following template installs node.js, python, pip, and other dependencies that are needed for the tool. 
More about definition file: [https://sylabs.io/guides/3.0/user-guide/definition_files.html](https://sylabs.io/guides/3.0/user-guide/definition_files.html)

6. Create .sif file: `sudo singularity build clinService.sif clinService.def` This process might take sometime. 
7. Move the newly created .sif file in `tool_bin` directory: `cp clinService.sif ../../tool_bin/clinService`
