This document outlines the process for adding a new backend service to gru.
In general, everyone is responsible for developing and testing their own
services, which will be submitted as pull requests and reviewed before
merging.

Here are the steps:

1. Set up a personal backend development environment on dev.backend.iobio.io by
   following the documentation [here](./gru_sandbox.md). This only needs to be
   done once. You can keep your environment updated and use it for future
   service development.

2. Propose your new service to Anders and/or Yi. Ultimately one of them will
   need to approve the pull request anyway, and this is a good opportunity to
   discuss architectural tradeoffs.

   It's very helpful if you bring the following information with you to this
   discussion:
   
   1. What language is your service written in? (prefer JavaScript when
      possible, using the [koajs framework](https://koajs.com/) (it's basically
      a more modern framework by the same people who made express).
   
   2. What do you expect a typical request to your service to look like? ie
      what are the inputs/outputs?

   3. What is the expected runtime for each request? For requests that take a
      long time to run, we might need to look into caching options or other
      mitigations.

   4. What data does your service depend on, and how often does it need to
      be updated? Small datasets are no problem. Large datasets that don't
      change often can be baked into our data volume. Large datasets that
      change often are a problem, and will need to be discussed.

3. Change your frontend code to point to your development backend instance.
   The URL will end up looking like `dev.backend.iobio.io:<PORT_NUM>`, ie
   `dev.backend.iobio.io:9004`, for example. Don't merge these frontend
   changes.  Frontend master branches should never point at
   `dev.backend.iobio.io`. 

   Note that most iobio functionality may run
   considerably slower on your development instance, because the production
   backend is load-balanced across multiple nodes.

4. Develop and test your backend on a git branch. JavaScript services can
   often be added directly to the codebase. Tools that require other runtimes
   usually need to be packaged into singularity containers. Ask Anders for
   details.

5. Submit a pull request.

6. Once the review process is done, Anders will merge your service and deploy
   it to production.
