README file for R package supporting the paper "scTSSR: gene expression recovery for single-cell RNA sequencing using two-side sparse self-representation".


The scTSSR package has the following R-package dependencies: SAVER, keras, tensorflow.
The dependent packages will be automatically installed along with scTSSR. You can use the following commands to install scTSSR from GitHub.

Installation
----------------------
### Step 1. If the devtools package has been not installed, install the devtools package first. 

`install.packages("devtools")`

### Step 2. Load the devtools package.

`library("devtools")`

### Step 3. Install the scTSSR package from GitHub.

`install_github("Zhangxf-ccnu/scTSSR")`


Useage
----------------------
Load the library scTSSR in R console, by running

`library(scTSSR)`

Taking the baron the dataset as an example, run the following code:

`data("baron")`

`baron\_imputation\_result <- scTSSR(baron$count.samp, percent=0.05, learning\_rate=0.0001, epochs=100)`

For detialed usages, please refer to "scTSSR-manual.pdf".

Codes for reproducing the three downstream analyses (such as differential expression analysis, cell clustering analysis and pseudotime analysis) are available in [scTSSR-scTSSR2_experiments_codes](https://github.com/Zhangxf-ccnu/scTSSR-scTSSR2_experiments_codes) and Zenodo website with [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423597.svg)](https://doi.org/10.5281/zenodo.6423597).


Contact
------------------------
Please do not hesitate to contact Miss Ke Jin (kej13@mails.ccnu.edu.cn) or Dr. Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.







