README file for R package supporting the paper "scTSSR: gene expression recovery for single-cell RNA sequencing using two-side sparse self-representation".


Contents of this archive
------------------------
This archive contains 

(1) pkg: subdirectory that contains the R package.

(2) scTSSR-manual.pdf: reference manual.

(3) trajectory_analysis: subdirectory that contains scripts for performing cell trajectory inference on the deng dataset (Section 3.5 in the main text). You can run the "trajectory analysis.R" file to carry out the experiments.

The scTSSR package has the following R-package dependencies: SAVER.
The dependent packages will be automatically installed along with scTSSR. You can use the following commands to install scTSSR from GitHub.


# Step 1. If the devtools package has been not installed, install the devtools package first. Invoke R and then type

install.packages("devtools")

# Step 2. Load the devtools package.

library("devtools")

# Step 3. Install the scTSSR package from GitHub.

install_github("Zhangxf-ccnu/scTSSR", subdir="pkg")


Useage

Load the library scTSSR in R console, by running

library(scTSSR)

Taking the baron the dataset as an example, run the following code:

data("baron") 

baron_imputation_result = scTSSR(baron$count.samp)

For detialed usages, please refer to "scTSSR-manual.pdf".

Please do not hesitate to contact Miss Ke Jin (kej13@mails.ccnu.edu.cn) or Dr. Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.
