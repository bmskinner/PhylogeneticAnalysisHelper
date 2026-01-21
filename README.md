# Introduction

This R package bundles common functions I have written for various phylogentic analyses around plotting and tree manipulation.

It also provides functions to invoke external binaries for standard phylogenetics packages like `muscle` and `iqtree`. It is assumed that all the binaries are already on the PATH; this package just uses R's `system` function to invoke them and return output file paths in a consistent manner.

# Installation

```
# Install devtools from CRAN if needed

devtools::install_githib("bmskinner/PhylogeneticsAnalysisHelper")
```
