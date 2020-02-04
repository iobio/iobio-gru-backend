This document will outline the steps to build a singularity container. 

**Building a singularity container requires the following:** 
* A Definition File 
* Singularity build command 




**Creating a Definition file:** 
A singularity definition file (.def) is a set of blueprints explaining how to build a custom container. Here are the steps to create a .def file.
1. Go to tools directory. example: `cd adit/iobio-gru-backend/tools` 
2. Create a directory for your tool. Let us assume that the name of our tool is clin service: `mkdir clinService`
3. Go to the newly created directory: `cd clinService`
4. Create a new .def file: `vi clinService.def`
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


Once a singularity container is built, create new routes so the front-end can make requests to your service. You will also need to create a .sh file so that your request can invoke the tool. 

**Adding routes** 

1. Go to src directory:  `cd adit/iobio-gru-backend/src`
2. Open index.js: `vi index.js` 
3.  Add new endpoints for your service. Example 
```javascript
// GET request 
router.get('/clinService', async (ctx) => {
  const args = [ctx.query.urlParamterName];
  //Add or update code.. 
  await handle(ctx, 'clinService.sh', args);
});
```  
You can look at other routes in the same file for the POST requests. 
