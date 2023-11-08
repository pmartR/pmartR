### November 14, 2023
- Skipping long-running tests to speed up checks.  These tests are still run with github actions.
- Update LICENSE.


### November 9, 2023
**Dealt with error relating to missing package `pmartRdata` in Suggests.**  Package moved to drat repository and DESCRIPTION updated to use `Additional_repositories`.  Examples now run conditionally on whether pmartRdata can be loaded with `requireNamespace` to prevent errors during checks.

**Added vignettes to .Rbuildignore to avoid warnings.**  Vignettes take a long time to run and are published to a github pages.

**Notes Fixed**
- Documentation lines with > 100 characters.
- Fixed documented but not used arguments.
- Removed C++11 specification
- Consolidated non-standard top level file `disclaimer.txt` into README
- Suppressed long running examples when running `--as-cran`.

**Remaining Notes**
In rhub checks, package size exceeds 5mb due to the 'lib' directory, which contains the compiled C++ code.  This is required for the package to function.  Previous submission did not appear to contain this note.

Note about pmartRdata not available for checks.  As mentioned this is because it is hosted in a drat repository.  This is only needed for vignettes and examples.

**Files being written to disk**
Fixed tests writing to the ~/Downloads folder.  It was also noted that there were other files being written to src/R/.svn/pristine/.... however we can't see how our code would write to this location, I'm unable to find these files on ubuntu/debian containers.
