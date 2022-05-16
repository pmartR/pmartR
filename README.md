# pmartR

[![DOI](https://zenodo.org/badge/69275428.svg)](https://zenodo.org/badge/latestdoi/69275428)

This R package provides functionality for quality control processing, statistical analysis and visualization of mass spectrometry (MS) omics data, in particular proteomic (either at the peptide or the protein level; isobaric labeled or unlabled), lipidomic, and metabolomic data. This includes data transformation, specification of groups that are to be compared against each other, filtering of feature and/or samples, data normalization, data summarization (correlation, PCA), and statistical comparisons of groups of interest (ANOVA and/or independence of missing data tests). Example data to be used with this packages can be found in [pmartRdata](https://github.com/pmartR/pmartRdata).


## Installation:

``` r
devtools::install_github("pmartR/pmartR")
```


## Tutorial:

To get started, see the package documentation and function reference located [here](https://pmartr.github.io/pmartR/index.html).

## Data:

Example peptide (both unlabeled and isobaric labeled), protein, metabolite and lipid data are available in the __pmartRdata__ package available on Github, [here](https://github.com/pmartR/pmartRdata)
 
## Citation:

To cite this package, please use the following:

Stratton KG, Webb-Robertson BJ, McCue LA, Stanfill B, Claborne D, Godinez I, Johansen T, Thompson AM, Burnum-Johnson KE, Waters KM, Bramer LM. pmartR: quality control and statistics for mass spectrometry-based biological data. Journal of proteome research. 2019 Jan 14;18(3):1418-25.

@article{stratton2019pmartr,
  title={pmartR: quality control and statistics for mass spectrometry-based biological data},
  author={Stratton, Kelly G and Webb-Robertson, Bobbie-Jo M and McCue, Lee Ann and Stanfill, Bryan and Claborne, Daniel and Godinez, Iobani and Johansen, Thomas and Thompson, Allison M and Burnum-Johnson, Kristin E and Waters, Katrina M and others},
  journal={Journal of proteome research},
  volume={18},
  number={3},
  pages={1418--1425},
  year={2019},
  publisher={ACS Publications}
}
