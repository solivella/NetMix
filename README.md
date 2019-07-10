[![Build Status](https://travis-ci.com/solivella/NetMix.svg?branch=master)](https://travis-ci.com/solivella/NetMix) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/NetMix)](https://cran.r-project.org/package=NetMix)

# NetMix
NetMix: Dynamic Mixed-Membership Network Regression Model. Variational EM estimation of mixed-membership stochastic blockmodel for networks, incorporating node-level predictors of mixed-membership vectors, as well as dyad-level predictors. For networks observed over time, the model defines a hidden Markov process that allows the effects of node-level predictors to evolve in discrete, historical periods. In addition, the package offers a variety of utilities for exploring results of estimation, including tools for conducting posterior predictive checks of goodness-of-fit and several plotting functions. 

# Installation
To install the latest CRAN version, copy and paste the following in a running R console:
`install.packages("NetMix")`.

For the latest development version, copy and paste the following (assuming you have the `devtools` package installed):
`devtools::install_github("solivella/NetMix")`
