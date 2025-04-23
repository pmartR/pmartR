# Contributing to pmartR

## Pull request process

### pmartR Contributors

If you have permissions to push branches directly to the pmartR repository, create your branch, make changes, and then open a pull request with a link to the issue it solves.  Specifically:

1. pull the latest changes from the master branch
```bash
git checkout master
git pull
```

2. create your branch
```bash
git checkout -b your-branch
```

3. make changes, run `devtools::test()` and/or `devtools::check()` and commit

4. push changes to the remote
```bash
git push -u origin your-branch
```

5. open a pull request in the github UI to the master branch.

**Make sure your branch can be fast-forwarded or cleanly rebased onto the master branch.**  In the case that this cannot be done automatically in the github UI, you can do this by rebasing locally and resolving any conflicts, then force-pushing your changes to your branch:

1. pull the latest changes from the master branch
```bash
git checkout master
git pull
```

2. rebase your branch onto the master branch
```bash
git checkout your-branch
git rebase -i master
```

The -i flag gives some control over how you want to rebase.  You can squash commits, edit commit messages, etc.  Omit the -i flag if you just want to do a default rebase attempt.  If you have conflicts, you will need to resolve them and then continue the rebase process.

3. run `devtools::test()` and/or `devtools::check()`

4. force push your branch to the remote
```bash
git push origin your-branch --force
```

We prefer using rebase to merge, as merge can create very messy, unintuitive history.  We are also assuming people are all branching off of master, so no one will have branched off of your branch before you rewrite the history with rebase.

### Outside pmartR team

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("pmartR/pmartR", fork = TRUE)`.

*   Install all development dependencies with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

## Testing and Deployment

### Local Checks

Whenever making significant code changes, run `devtools::test()` and ensure that all tests pass.  A more complete check of the code can be run using `devtools::check()`.

### Running Automatic Checks

We use automatic checks via github actions to make sure code is formatted correctly, tests pass, examples run, etc.  These run when a pull request to master is created/edited or a manual run is triggered.  See the documentation for [github actions](https://docs.github.com/en/actions) for more info and the `.github/workflows` folder for the yaml files that specify the actions to be performed and when they will be run.  Manual triggers of workflows can be done as described [here](https://docs.github.com/en/actions/managing-workflow-runs/manually-running-a-workflow)

We include a github actions workflow to optionally run `rhub` (https://github.com/r-hub/rhub) checks.  Depending on where the branch you want to test is, you can run:

```r
rhub::rhub_check(gh_url = 'https://github.com/pmartR/pmartR', branch = 'your-branch')
```

Where you can replace the `gh_url` with the address of your fork if you are working from one.  This will prompt you to, among other things, select several environments/containers (e.g. ubuntu, arm64 mac) to run CRAN-like checks on.  Once completed, a series of jobs will spin up under the 'Actions' tab of the pmartR github repo or your fork if you ran it there.

### pkgdown deploy site
Pull requests and pushes to branches specified in .github/workflows/pkgdown.yaml will trigger a build of a web page containing function documentation and rendered vignettes.  This requires obviously that your vignettes build without error, and that your documentation is valid.  The site will **deploy** only when we do a release on github.

## Submitting to CRAN

After checks have passed (preferably the `rhub` checks), the easiest way to submit the package to CRAN for review is using `devtools`.  First, update `cran-comments.md` specifying any changes regarding CRAN checks (e.g. fixed NOTE/WARNINGs or other things that they would care about when deciding to accept/reject a submission).  Then we call `devtools::submit_cran(args = '--no-build-vignettes')`.  You will be prompted with a series of questions and then the package will be built and submitted.  The comments in `cran-comments.md` will be automatically included in the submission.

The package maintainer as listed in DESCRIPTION will be sent an email confirming the submission.  Once they confirm it will be sent to the CRAN maintainers for review.

Alternatively, you can use `R CMD build --no-build-vignettes pmartR` to build the package tarball and then manually upload it to the CRAN [submission site](https://cran.r-project.org/submit.html).

## Documentation and tests

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  
*  Vignettes should be put under `/vignettes` and should be written in `Rmarkdown`.
*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  
