# Estimating Network-Mediated Causal Effects via Principal Components Network Regression

We develop a method to decompose causal effects on a social network into an indirect effect mediated by the network, and a direct effect independent of the social network. To handle the complexity of network structures, we assume that latent social groups act as causal mediators. We develop principal components network regression models to differentiate the social effect from the non-social effect. Fitting the regression models is as simple as principal components analysis followed by ordinary least squares estimation. We prove asymptotic theory for regression coefficients from this procedure and show that it is widely applicable, allowing for a variety of distributions on the regression errors and network edges. We carefully characterize the counterfactual assumptions necessary to use the regression models for causal inference, and show that current approaches to causal network regression may result in over-control bias. The method is very general, so that it is applicable to many types of structured data beyond social networks, such as text, areal data, psychometrics, images and omics.

## To replicate our computational results

We use [`renv`](https://rstudio.github.io/renv/) to record package dependencies and [`targets`](https://books.ropensci.org/targets/) to coordinate our simulation study and data analysis.

To replicate our results, clone this Github repository to your local computer. Once you have the repository cloned locally, re-create the project library by calling

``` r
# install.packages("renv")
renv::restore()
```

Note that this will install a development version of [`latentnetmediate`](https://github.com/alexpghayes/latentnetmediate). `latentnetmediate` implements the principal components network regression methods, and, as of the time of this writing, is functional but undocumented. `latentnetmediate` also contains several datasets used in the paper.

At this point, you should be ready to replicate our simulation and performance comparison results. The computational work is organized into distinct `targets` projects (see [here](https://books.ropensci.org/targets/projects.html#multiple-projects)) for details. You will need to build these one at a time:

``` r
Sys.setenv(TAR_PROJECT = "addhealth")
tar_make()

Sys.setenv(TAR_PROJECT = "glasgow")
tar_make()

Sys.setenv(TAR_PROJECT = "healthyminds")
tar_make()

Sys.setenv(TAR_PROJECT = "misspecification")
tar_make()

Sys.setenv(TAR_PROJECT = "simulations")
tar_make()
```

This is a computationally intensive project. It's likely that, at some point, the targets pipeline will crash for some sundry computational reason. If this is the case, simply re-run `tar_make()`. `tar_make()` will only attempt to re-run incomplete portions of the build pipeline.

To fully replicate the paper, increase the number of replications in the misspecification and simulation studies by setting:

- `tar_target(num_chunks, 1)` -> `tar_target(num_chunks, 10)` on line 73 of `_misspecification.R` and
- `tar_target(num_chunks, 1, deployment = "main")` -> `tar_target(num_chunks, 30, deployment = "main")` on line 196 of `_simulations.R`

otherwise the misspecification and simulation studies will run using fewer replicates than in the paper (for computational purposes). We recommend you make sure the code runs with the default values before increasing them to the full replication size.
