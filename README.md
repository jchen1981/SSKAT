# SSKAT
<ins>S<ins>mall-<ins>S<ins>ample <ins>K<ins>ernel Machine <ins>A<ins>ssociation <ins>T<ins>est v1.0

## Overview
The classic SKAT (SNP-set/sequence kernel association test) could lose power if (1) the sample size is small, or (2) the coefficient of variation of the kernel spectrum is small. The power loss is due to ignoring the uncertainty in the error variance estimate.  We propose several versions of adjusted SKAT for different types of outcomes to improve the power of SKAT for its unfavorable situations.  The adjusted SKAT circumvents the difficulty of estimating the error variance (scale) under small sample sizes by deriving a scale-free statistic similar to F-statistic.  We implement the small-sample adjusted SKAT for

* Univariate continuous and binary outcome (aSKAT)
* Multivariate continuous outcome (mSKAT)
* Correlated continuous outcome (cSKAT)
* Generalized mixed effects model (VSAT)


## Reference

* Chen J, et al. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome Association Studies. Genet Epidemiol. 40(1):5-19.
* Zhan X, et al. (2017)  A small-sample multivariate kernel machine test for microbiome association studies. Genet Epidemiol. 41(3):210-220.
* Zhan X, et al. (2018) A small-sample kernel association test for correlated data with application to microbiome association studies. Genet Epidemiol. 42(8):772-782. 
* Zhan X, et al. (2020) Variant-Set Association Test for Generalized Linear Mixed Model. Genet Epidemiol. Submitted



## Installation         

```
install.packages(c("Matrix", "lme4qtl", "CompQuadForm", "lme4"))
install.packages("devtools")
devtools::install_github("variani/lme4qtl")
devtools::install_github("jchen1981/SSKAT")
```



## Examples

### Univariate continuous and binary outcome
```
library("SSKAT")
set.seed(123)
Y <- rnorm(100)
Z <- matrix(rnorm(200), 100, 2)
G <- matrix(rbinom(1000, 1, 0.25), 100, 10)
K <- G %*% t(G)
data <- data.frame(Y = Y, Z = Z)
aSKAT(Y ~ Z, data = data, K = K, type = 'continuous')
aSKAT(Y ~ 1, K = K, type = 'continuous')  # No covariate
Y <- rbinom(100, 1, 0.5)
data <- data.frame(Y = Y, Z = Z)
aSKAT(Y ~ Z, data = data, K = K, type = 'binary')
aSKAT(Y ~ 1, data = data, K = K, type = 'binary')    # No covariate
```


### Multivariate continuous outcome with de-correlation
```
library("SSKAT")
set.seed(123)
Y <- matrix(rnorm(400), 100, 4)
Z <- matrix(rnorm(200), 100, 2)
G <- matrix(rbinom(1000, 1, 0.25), 100, 10)
K <- G %*% t(G)
data <- data.frame(Z = Z)
mSKAT(~ Z, data = data, Y = Y, K = K) 
mSKAT(Y = Y, K = K)    # No covariate
mSKAT(~ 1, Y = Y, K = K)    # No covariate
```

### Multivariate continuous outcome without de-correlation
```
library("SSKAT")
set.seed(123)
L <- matrix(rnorm(1000), 100, 10)  # Latent factor to induce correlation in Y
Y <- scale(L %*% matrix(rnorm(10 * 40), 10, 40)) + matrix(rnorm(100 * 40), 100, 40)
G <- scale(Y %*% matrix(rnorm(40 * 10), 40, 10)) * 0.2 + matrix(rnorm(100 * 10), 100, 10)
K <- G %*% t(G)
mSKAT2(Y = Y, K = K)  # No de-correlation 
```


### Univariate correlated continuous outcome
```
library("SSKAT")
set.seed(123)
Y <- rnorm(100)
Z <- matrix(rnorm(200), 100, 2)
ID <- gl(20, 5)
G <- matrix(rbinom(1000, 1, 0.05), 100, 10)
K <- G %*% t(G)
data <- data.frame(Y = Y, Z = Z, ID = ID)
cSKAT(Y ~ Z + (1 | ID), data = data, K = K)
cSKAT(Y ~ (1 | ID), data = data, K = K)  # No covariate
```

### Univariate correlated general outcome (Generalized Mixed Effects Model, GLMM)
```
library("SSKAT")
set.seed(123)
m <- 20; k <- 4; n <- m * k; s <- 20
x <- rnorm(n); z <- matrix(sample(0 : 1, n * s, repl = TRUE, prob = c(0.8, 0.2)), n, s)
a <- rnorm(s); b <- 1
y <- rbinom(n, 1, prob = binomial()$linkinv(z %*% a + x * b + 0.5 * rnorm(m)[gl(m, k)]))
data <- data.frame(y = as.numeric(y), x = x, subject = gl(m, k))
rownames(data) <- paste(1:n)
K <- z %*% t(z)
VSAT(formula.H0 = y ~ x + (1 | subject), data = data, K = K, family = binomial()) 
VSAT(formula.H0 = y ~ (1 | subject), data = data, K = K, family = binomial()) # No covariate
```

