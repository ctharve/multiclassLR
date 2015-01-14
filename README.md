multiclassLR: multinomial regression in R  
============
This package contains an implementation of multinomial logistic regression and is optimized with the optim function (this is quite slow, so alternative methods are currently being developed). It can be used for any data set with both categorical and continuous variables, but is specifically designed to be used on the output of MESH ([Wen et al, 2012]). For a set of SNPs, each of which has been subjected to a treatment and corresponding control, MESH estimates posterior probabilites of four configurations: treatment only Allele Specific Expression (ASE), control only ASE, ASE in both treatment/control, & ASE in neither treatment/control. These posteriors comprise the 4 classes (labels) to be used in the multiclass regression.
 
Each individual SNP X treatment combination is labeled by the configuration with maximal posterior probability. In terms of features to regress over we consider the abundance of transcripts overlapping each SNP, indicators for the treatment, and indicators for the cell-line used in the experiment. 

## 1. Installation

To install from within an R session:

```R
require(devtools)
install_github('ctharve/multiclassLR')
library('multiclassLR')
```

## 2. Data Preparation

<!-- links -->
[Wen et al, 2012]:http://arxiv.org/abs/1111.1210
