# pmartRqc

This R package provides functionality for quality control processing of mass spectrometry (MS) omics data, in particular proteomic (either at the peptide or the protein level), lipidomic, and metabolomic data. This includes data transformation, specification of groups that are to be compared against each other, filtering of feature and/or samples, data normalization, and data summarization (correlation, PCA). Example data to be used with this packages can be found in [pmartRdata](https://github.com/pmartR/pmartRdata).


``` r
devtools::install_github("pmartR/pmartR")
```

For Mac Users, to enable the C++11 plugin, run the following in your R console before installing this package:
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
 
