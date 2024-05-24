# pmartR
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/69275428.svg)](https://zenodo.org/badge/latestdoi/69275428)
[![R-CMD-check](https://github.com/pmartR/pmartR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pmartR/pmartR/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/pmartR)](https://CRAN.R-project.org/package=pmartR)
<!-- badges: end -->

This R package provides functionality for quality control processing, statistical analysis and visualization of mass spectrometry (MS) omics data, in particular proteomic (either at the peptide or the protein level; isobaric labeled or unlabled), lipidomic, and metabolomic data. This includes data transformation, specification of groups that are to be compared against each other, filtering of feature and/or samples, data normalization, data summarization (correlation, PCA), and statistical comparisons of groups of interest (ANOVA and/or independence of missing data tests). Example data to be used with this packages can be found in [pmartRdata](https://github.com/pmartR/pmartRdata).


## Installation:

This package makes use of several packages hosted on BioConductor.  If you are encountering warnings about unavailable BioConductor packages such as `pcaMethods`, you may need to add them to `options("repos")`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

options("repos" = BiocManager::repositories())
```

(**Recommended**) Install from CRAN:
```r
install.packages("pmartR")

# or 

BiocManager::install("pmartR")
```

If you are on Mac/Windows and have a recent R version, you can skip compilation by installing from binaries, see the [pmartR CRAN page](https://cran.r-project.org/web/packages/pmartR/index.html) for available binaries.

```r
install.packages("pmartR", type = "binary")

# or 

BiocManager::install("pmartR", type = "binary")
```


To install the latest *release*:
```r
devtools::install_github("pmartR/pmartR@*release")
```

To install a specific release, say v2.4.0:

```r
devtools::install_github("pmartR/pmartR@v2.4.0")
```

**(Not recommended, since these changes are likely still being tested)** You can also install the latest changes to master:

```r
devtools::install_github("pmartR/pmartR")
```

### Problems with rcppArmadillo and gfortran on mac

There is a problem that causes pmartR to fail compiling cpp code, which has something to do with rcppArmadillo and certain installations of gfortran.  See these posts that try to explain the issue:  [1](https://stackoverflow.com/questions/64992467/mac-clang-installation-seems-to-override-gcc-install) [2](https://stackoverflow.com/questions/29992066/rcpp-warning-directory-not-found-for-option-l-usr-local-cellar-gfortran-4-8/29993906#29993906) [3](https://forum.posit.co/t/setting-up-travis-ci-on-linux-with-an-r-package-that-uses-rcpparmadillo/53910/3).  The simplest solution if you are on Mac/Windows and have a recent R version is to install from pre-built binaries (see installation section).  If you cannot install from binaries, two solutions we have found are:

1.  Install gfortran from a recommended source (not homebrew): 
    - [This CRAN-approved resource for build tools on mac](https://mac.r-project.org/tools/) lists two versions of gfortran and how to install them.
    - On Catalina 10.15.7 I downloaded and installed gfortran 8.2 from the link provided in [this blog post](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/#google_vignette)  
2.  When using the homebrew gfortran installation, add the line **FLIBS = -L\`gfortran -print-file-name=libgfortran.dylib | xargs dirname\`** to ~/.R/Makevars (a plain text file with no extention)

### gfortran and Apple silicon (M1/M2 chips)

There are similarly issues with compilation in newer Mac chips.  We recommend to install gcc-13 from homebrew `brew install gcc` or the universal version from https://mac.r-project.org/tools/.  

Additionally, some users experience errors with `ld:  Assertion failed ...` as seen [here](https://developer.apple.com/forums/thread/737707).  One solution is to use the old linker by making sure gcc uses the flag `-ld64` [(Xcode docs)](https://developer.apple.com/documentation/xcode-release-notes/xcode-15-release-notes#Linking).  To do this, you can edit `~/.R/Makevars` to include this flag, for example by appending it to LDFLAGS with +=:

```
# in ~/.R/Makevars
LDFLAGS+=-ld64
```

or specifying it in your compiler command:

```
# in ~/.R/Makevars
CC=/usr/local/bin/gcc -ld64
```

## Tutorial:

To get started, see the package documentation and function reference located [here](https://pmartr.github.io/pmartR/).

## Data:

Example peptide (both unlabeled and isobaric labeled), protein, metabolite and lipid data are available in the __pmartRdata__ package available on Github, [here](https://github.com/pmartR/pmartRdata)

## Contributing
See the [contributing docs](.github/CONTRIBUTING.md).

## Citation:

To cite this package, please the following:

Degnan, D. J.; Stratton, K. G.; Richardson, R.; Claborne, D.; Martin, E. A.; Johnson, N. A.; Leach, D.; Webb-Robertson, B.-J. M.; Bramer, L. M. PmartR 2.0: A Quality Control, Visualization, and Statistics Pipeline for Multiple Omics Datatypes. J. Proteome Res. 2023, 22 (2), 570–576. https://doi.org/10.1021/acs.jproteome.2c00610.

**BibTex:**

```
@article{degnan2023pmartr,
  title={pmartR 2.0: A Quality Control, Visualization, and Statistics Pipeline for Multiple Omics Datatypes},
  author={Degnan, David J and Stratton, Kelly G and Richardson, Rachel and Claborne, Daniel and Martin, Evan A and Johnson, Nathan A and Leach, Damon and Webb-Robertson, Bobbie-Jo M and Bramer, Lisa M},
  doi = {10.1021/acs.jproteome.2c00610},
  journal={Journal of Proteome Research},
  year={2023},
  publisher={ACS Publications}
}
```

## Disclaimer:

This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

      PACIFIC NORTHWEST NATIONAL LABORATORY
      operated by BATTELLE for the
      UNITED STATES DEPARTMENT OF ENERGY
      under Contract DE-AC05-76RL01830
