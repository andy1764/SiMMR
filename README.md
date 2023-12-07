# Similarity-based multimodal regression (SiMMR)

**Maintainer**: Andrew Chen, chenandr@musc.edu

**License**: Artistic License 2.0

SiMMR is a distance-based approach for analysis of multiple data modalities. Several test statistics are proposed, which are all implemented in the `simmr` function. Further details are available in our preprint (see [Citations](#3.-citations)).
 
## 1. Installation
The R package can be installed via devtools by running the following code

```
# install.packages("devtools")
devtools::install_github("andy1764/SiMMR", build_vignettes = TRUE)
```

## 2. Usage
SiMMR can be run by first computing distance matrices for each modality, then inputting as a list of distance matrices to `simmr`. Below is an example call using simulated data:

```
D <- list(dist(rnorm(10)), dist(rnorm(10)), dist(rnorm(10)))
X <- list("var" = runif(10))
simmr(D, X, "var")
```

## 3. Citations
Please cite the following article:

> Chen, A. A., Weinstein, S. M., Adebimpe, A., Gur, R. C., Gur, R. E., Merikangas, K. R., Satterthwaite, T. D., Shinohara, R. T., & Shou, H. (2023). Similarity-based multimodal regression. *Biostatistics*, kxad033. https://doi.org/10.1093/biostatistics/kxad033

