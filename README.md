# pmartR

This R package provides functionality for quality control processing, statistical analysis and visualization of mass spectrometry (MS) omics data, in particular proteomic (either at the peptide or the protein level), lipidomic, and metabolomic data. This includes data transformation, specification of groups that are to be compared against each other, filtering of feature and/or samples, data normalization, and data summarization (correlation, PCA). Example data to be used with this packages can be found in [pmartRdata](https://github.com/pmartR/pmartRdata).


## Installation:

``` r
devtools::install_github("pmartR/pmartR")
```

It’s been noted, here https://stackoverflow.com/questions/51257009/is-rtools-incompatible-with-r-version-3-5-1 , that there are issues with Rtools for R version 3.5.1.
 
If you run into the same problem when trying to install using `devtools::install_github()` and the suggested fix in the above link does not work, you can clone or download the pmartR library to your computer and install it locally.

It's also been noted, here https://github.com/rcppsmc/rcppsmc/issues/27 , that there are issues with xcode. This manifests itself as an error stating "math.h" not found when trying to install and restart a package. The following fix, at the terminal, has worked for us:

xcode-select --install


## Tutorial:

To get started, see the package documentation and function reference located [here](https://pmartr.github.io/pmartR/index.html).

## Data:

Example peptide, protein, metabolite and lipid data are available in the __pmartRdata__ package available on Github, [here](https://github.com/pmartR/pmartRdata)
 
