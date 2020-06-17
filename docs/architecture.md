# High-level components

* gru (this repo)
* Phenolyzer
* fibridge (local file proxy)
* Multi-align 
* stealthcheck


## gru

This is the primary component of the iobio backend. Currently, gru has 2 jobs:

1. Remote invocation of command line bioinformatics pipelines, which combines
   tools such as samtools, VEP, bcftools, etc.
2. Several other services, including /geneinfo, /gene2pheno, and /hpo.
3. Serving static files such as ClinVar data

Eventually everything other than the pipelines will likely be moved into
iobio-backend-services alongside phenolyzer (see below).


## Phenolyzer

Phenolyzer is currently the (only) subcomponent of [iobio-backend-services][0]. It
is a basic services. It takes a phenotype `term` query parameter, and returns a
list of genes associated with the term.

The phenolyzer services maintains a simple flat-file cache of previously
requested terms, in `./cache`.


## fibridge

[fibridge][3] is a generic service which provides a way to proxy HTTP
connections to a WebSocket source. iobio uses it for local file support. The
way it works is the iobio app opens a file on the user's machine (for example a
bam file) and also opens a WebSocket channel to the fibridge server. The
fibridge server provides a URL for the file. The app then invokes the gru
backend, giving it the fibridge URL.  From gru's perspective, it's just an HTTP
request, but in reality the app is on the other end, serving the file over the
WebSocket, proxied through the fibridge server.


## Multialign

[multialign][1] is a separate service for doing multi-way alignments. It will likely
eventually be integrated into iobio-backend-services.


## stealthcheck

[stealthcheck][2] is a simple health check/restart/email alert service.
Currently it keeps an eye on fibridge, phenolyzer, and a few other things which
are not critical for running the iobio backend.


# gru vs other services

It may be confusing why some functionality is integrated into gru, and other
functionality isn't. This section attempts to explain this.

This first thing to note is that gru has evolved over time, so some things have
been added and others will likely be removed/moved elsewhere.

Other than that, there are a coupld major concerns for the iobio backend:

* Security (ie HIPAA compliance)
* Nature of different requests (ie cacheable vs not, request/response vs RPC)


## Security

This is the primary distinction. Some functionality provided by the iobio
backend is public in nature. For example, converting a phenotype into a list
of genes (phenolyzer), or serving ClinVar data. Other than performance
concerns, it doesn't really matter where these services are hosted.

However, any service which touches user-provided data must be able to run in
an environment that offers the same level of security and privacy that the
data itself has. This precludes for example running such tasks on AWS Lambda.

gru's job is to handle all of these types of requests.


## Request type

Most of the pipeline requests are really remote procedure calls (RPCs). They
even use/abuse HTTP POST requests with content type text/plain (event though
they're sending JSON) as a hack to get around the limitations of HTTP/CORS.
These requests are not cacheable. Many of the requests are also
resource-intensive, therefore gru is designed to be load-balanced across as
many nodes as necessary.

Many of the other requests are just generic HTTP requests, and are cacheable.
These are a much better fit for running on a single server, and relying on
caching for scalability. /gene2pheno, /geneinfo, phenolyzer, and static data
all fall into this category. The current plan is to migrate all these services
into a single cacheable service.


[0]: https://github.com/iobio/iobio-backend-services

[1]: https://github.com/iobio/multialign 

[2]: https://github.com/anderspitman/stealthcheck

[3]: https://github.com/anderspitman/fibridge-proxy-rs
