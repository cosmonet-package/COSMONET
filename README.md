# COSMONET package
COSMONET is a new R package that implements a novel multistage computational-statistical procedure combining screening techniques and network methods. It is able to identify prognostic gene signatures and predict patient survival outcome [1,2]. We applied (i) screening approaches to reduce the initial dimension from an high-dimensional space p to a moderate scale d and (ii) network-penalized Cox-regression approaches to model observed survival time by using genome-wide omic profiles while accounting for coordinated functioning of genes in the form of biological pathways or networks. We illustrate the use of our package by applying it to breast cancer datasets downloaded from TCGA portal. 

## Installation

### On Mac OS X

Open terminal and run ````xcode-select --install```` to install the command line developer tools.

### Installing from github

To install COSMONET via GitHub enter the following commands: 

````
install.packages("devtools")
devtools::install_github("cosmonet-package/COSMONET", repos=BiocManager::repositories())
````

Note. Execute the following steps:

1. Open Terminal
2. Write or paste in: ````defaults write org.R-project.R force.LANG en_US.UTF-8````
3. Close Terminal (including any RStudio window)
4. Start R

if R gives the below message when you install the R packages:

````
Error: (converted from warning) Setting LC_CTYPE failed, using "C"
Execution halted
Error in ...
````

### To get started

Once installed, upload COSMONET package with the following command:

````
library(COSMONET)
````

For more details, please consult vignette documentation.

## References

[1] Iuliano et al. (2018). Combining pathway identification and breast cancer survival prediction via screening-network methods. Frontiers in genetics, 9.

[2] Iuliano et al. (2016). Cancer markers selection using network-based Cox regression: A methodological and computational practice. Frontiers in physiology, 7, 208.
