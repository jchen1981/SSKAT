# Title: Small-sample Kernel Machine Association Test
# Version: 0.1
# Author: Jun Chen (chen.jun2@mayo.edu)
# Date: 2018/02/07 - 2020/08/01
sqrt.inv <- function (V2) {
	eig.obj <- eigen(V2, symmetric = TRUE)
	vectors <- eig.obj$vectors
	values <- eig.obj$values
	ind <- values >= 1e-10
	values <- values[ind]
	vectors <- vectors[, ind]
	
	temp <- t(vectors) / (values)
	Vi2 <- vectors  %*% temp
	
	temp <- t(vectors) / sqrt(values)
	Vi <- vectors  %*% temp
	
	return(list(Vi = Vi, Vi2 = Vi2, rank = length(values)))
}


### The following three functions are copied directly from 'MSKAT' package by Baolin Wu
### http://www.biostat.umn.edu/~baolin/research/mskat_Rcode.html
saddle = function(x,lambda){
	d = max(lambda)
	lambda = lambda/d
	x = x/d
	k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
	kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
	kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
	n = length(lambda)
	if (any(lambda < 0)) {
		lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
	} else if (x>sum(lambda)){
		lmin = -0.01
	} else {
		lmin = -length(lambda)/(2*x)
	}
	lmax = min(1/(2*lambda[lambda>0])) * 0.99999
	hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
	w = sign(hatzeta)*sqrt(2*(hatzeta * x-k0(hatzeta)))
	v = hatzeta*sqrt(kpprime0(hatzeta))
	if(abs(hatzeta)<1e-4){
		return(NA)
	} else{
		return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
	}
}

Liu.pval = function(Q, lambda){
	c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
	muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
	s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
	if(s1^2 > s2){
		a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
	} else {
		l = 1/s2;  a = sqrt(l);  d = 0
	}
	muX = l+d;  sigmaX = sqrt(2)*a
	
	Q.Norm = (Q - muQ)/sigmaQ
	Q.Norm1 = Q.Norm*sigmaX + muX
	pchisq(Q.Norm1, df = l,ncp=d, lower.tail=FALSE)
}

Sadd.pval = function(Q.all,lambda){
	sad = rep(1,length(Q.all))
	if(sum(Q.all>0)>0){
		sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
	}
	id = which(is.na(sad))
	if(length(id)>0){
		sad[id] = Liu.pval(Q.all[id], lambda)
	}
	return(sad)
}

#Compute the tail probability of 1-DF chi-square mixtures
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
	pval = rep(0, length(Q.all))
	i1 = which(is.finite(Q.all))
	for(i in i1){
		tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
		if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
	}
	return(pval)
}
################################################################################################


aSKAT.c <- function (formula.H0, data = NULL, K, acc = 0.00001, lim = 10000, tol = 1e-10) {
	
	res <- resid(lm(formula.H0, data))
	s2 <- sum(res^2) 
	X1 <- model.matrix(formula.H0, data)
	
	D0  <- diag(length(res))  
	P0  <- D0 - X1 %*% solve(t(X1) %*% X1) %*% t(X1)
	PKP <- P0 %*% K %*% P0
	
	q <- as.numeric(res %*% K %*% res / s2)
	
	ee <- eigen(PKP - q * P0, symmetric = T) 
	
	lambda <- ee$values[abs(ee$values) >= tol]

	p.value <- KAT.pval(0, lambda=sort(lambda, decreasing=T), acc = acc, lim = lim)
	
	return(list(p.value=p.value, Q.adj=q))
}


aSKAT.b <- function (formula.H0, data = NULL, K, acc = 0.00001, lim = 10000, tol = 1e-10) {
	
	X1 <- model.matrix(formula.H0, data)
	lhs <- formula.H0[[2]]
	y <- eval(lhs, data)
	
	y <- factor(y)
	

	if (nlevels(y) != 2) {
		stop('The phenotype is not binary!\n')
	} else {
		y <- as.numeric(y) - 1
	}
	
	glmfit <- glm(y ~ X1 - 1, family = binomial)
	
	betas <- glmfit$coef
	mu  <- glmfit$fitted.values
	eta <- glmfit$linear.predictors
	res.wk <- glmfit$residuals
	res <- y - mu
	
	w   <- mu * (1-mu)
	sqrtw <- sqrt(w)
	
	adj <- sum((sqrtw * res.wk)^2) 
	
	DX12 <- sqrtw * X1
	
	qrX <- qr(DX12)
	Q <- qr.Q(qrX)
	Q <- Q[, 1:qrX$rank, drop=FALSE]
	
	P0 <- diag(nrow(X1)) - Q %*% t(Q)
	
	DKD <- tcrossprod(sqrtw) * K
	tQK <- t(Q) %*% DKD
	QtQK <- Q %*% tQK 
	PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
	q <- as.numeric(res %*% K %*% res) / adj
	ee <- eigen(PKP - q * P0, symmetric = T, only.values=T)  		
	lambda <- ee$values[abs(ee$values) >= tol]
	
	p.value <- KAT.pval(0, lambda=sort(lambda, decreasing = T), acc = acc, lim = lim) 
	
	return(list(p.value=p.value, Q.adj = q))
}


#' Small-sample adjusted SKAT for a univariate binary/Gaussian outcome ('a' stands for 'adjusted') 
#'
#' Compute the adjusted score statistic and p-value for a univariate bianry/Gaussian outcome.
#' @param formula.H0  a two-sided linear formula object under the null, indicating the variables to adjust.
#' @param data  a data frame (required) containing the variables named in formula.
#' @param K    the kernel matrix, which quantifies the similarities between samples.
#' @param type   the outcome type: binary or continous. 
#' @param lim maximum number of integration terms for \code{davies} method. 
#' @param acc error bound for \code{davies} method. 
#' @param tol the eigenvalue cutoff, below which is considered to be 0. This is used to reduce the computation burden.
#' @return
#' \describe{A list containing
#'   \item{p.value}{ association p-value}
#'   \item{Q.adj}{ adjusted score statistic}
#' }
#' @keywords SKAT
#' @importFrom stats as.formula  model.matrix lm resid uniroot pnorm rnorm rbinom binomial glm pchisq
#' @import CompQuadForm
#' @rdname aSKAT
#' @author Jun Chen
#' @export
#' @references
#' Chen et al. (2016) Small Sample Kernel Association Tests for Human Genetic and Microbiome Association Studies.
#' Genet Epidemiol. 40(1):5-19.
#' @examples
#' set.seed(123)
#' Y <- rnorm(100)
#' Z <- matrix(rnorm(200), 100, 2)
#' G <- matrix(rbinom(1000, 1, 0.25), 100, 10)
#' K <- G %*% t(G)
#' data <- data.frame(Y = Y, Z = Z)
#' aSKAT(Y ~ Z, data = data, K = K, type = 'continuous')
#' aSKAT(Y ~ 1, K = K, type = 'continuous')  # No covariate
#' Y <- rbinom(100, 1, 0.5)
#' data <- data.frame(Y = Y, Z = Z)
#' aSKAT(Y ~ Z, data = data, K = K, type = 'binary')
#' aSKAT(Y ~ 1, data = data, K = K, type = 'binary')    # No covariate


aSKAT <- function(formula.H0, data = NULL, K, type = c('binary', 'continuous'), acc = 0.00001, lim = 10000, tol = 1e-10) {
	type <- match.arg(type)
	
	if (type == 'continuous') {
		res <- aSKAT.c(formula.H0, data, K, acc = acc, lim = lim, tol = tol) 
	}
	
	if (type == 'binary') {
		res <- aSKAT.b(formula.H0, data, K, acc = acc, lim = lim, tol = tol) 
	}
	
	res
}


#' Small-sample adjusted SKAT for a multivariate (continuous) outcome ('m' stands for 'multivariate') with de-correlation. 
#' De-correlation is effective when the number of outcome variables is small compared to the sample size.
#' 
#' Compute the adjusted score statistic and p-value for a multivariate continuous outcome with de-correlation.
#' @param formula.H0  a one-sided linear formula object under the null, indicating the variables to adjust.
#' @param data  a data frame containing the variables named in formula. It is required if the formula is not null or the intercept model (~ 1).
#' @param Y    a matrix for the multivariate outcome, row - samples, col - phenotypes.
#' @param K    the kernel matrix, which quantifies the similarities between samples.
#' @param lim maximum number of integration terms for \code{davies} method. 
#' @param acc error bound for \code{davies} method. 
#' @param tol the eigenvalue cutoff, below which is considered to be 0. This is used to reduce the computation burden.
#' 
#' @return A list containing
#' \describe{
#'   \item{p.value}{ association p-value}
#'   \item{Q.adj}{ adjusted score statistic}
#' }
#' @keywords SKAT
#' @importFrom stats as.formula  model.matrix lm resid uniroot pnorm rnorm rbinom binomial cov
#' @import CompQuadForm
#' @rdname mSKAT
#' @author Jun Chen
#' @export
#' @references
#' Zhan X, et al. (2017)  A small-sample multivariate kernel machine test for microbiome association studies.
#' Genet Epidemiol. 41(3):210-220.
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(400), 100, 4)
#' Z <- matrix(rnorm(200), 100, 2)
#' G <- matrix(rbinom(1000, 1, 0.25), 100, 10)
#' K <- G %*% t(G)
#' data <- data.frame(Z = Z)
#' mSKAT(~ Z, data = data, Y = Y, K = K) 
#' mSKAT(Y = Y, K = K)    # No covariate
#' mSKAT(~ 1, Y = Y, K = K)    # No covariate

mSKAT <- function(formula.H0 = NULL, data = NULL, Y, K, acc = 0.00001, lim = 10000, tol = 1e-10) {
	
	n <- ncol(K)
	if (!is.null(formula.H0)) {
		X <- model.matrix(formula.H0, data = data)
		if (nrow(X) == 0) {
			X <- rep(1, n)
		}
	} else {
		X <- rep(1, n)
	}
	
	Y <- resid(lm(Y ~ X - 1))
	V2 <- cov(Y)
	
	obj <- sqrt.inv(V2)
	Vi <- obj$Vi
	p <- obj$rank
	
	Y0 <- Y %*% Vi
	
	I <- diag(n)
	P <- I - X %*% solve(t(X) %*% X) %*% t(X)
	
	t0 <-  sum(Y0 %*% t(Y0) * K) / sum(Y0 * Y0)
	
	lambda1 <- eigen(P %*% K %*% P - t0 * P, only.values = T)$values
	lambda1 <- lambda1[abs(lambda1) >= tol]
	
	lambda2 <- rep(1, p)
	lambdas <- as.vector(outer(lambda2, lambda1, '*'))
	
#	davies(0, lambda=sort(lambdas, decreasing=T))$Qq
	p.value <- KAT.pval(0, lambda=sort(lambdas, decreasing=T), acc = acc, lim = lim) 
	list(p.value = p.value, Q.adj = t0)
}


#' Small-sample adjusted SKAT for a multivariate (continuous) outcome ('m' stands for 'multivariate') without de-correlation. 
#' De-correlation could be problematic when the number of outcome variables is not small compared to the sample size,
#' due to the large estimation error for the covariance matrix. In such case, de-correlation is not recommended unless
#' the correlation is very high. 
#' 
#' Compute the adjusted score statistic and p-value for a multivariate continuous outcome without de-correlation.
#' @param formula.H0  a one-sided linear formula object under the null, indicating the variables to adjust.
#' @param data  a data frame containing the variables named in formula. It is required if the formula is not null or the intercept model (~ 1).
#' @param Y    a matrix for the multivariate phenotypes, row - samples, col - phenotypes.
#' @param K    the kernel matrix, which quantifies the similarities between samples.
#' @param lim maximum number of integration terms for \code{davies} method. 
#' @param acc error bound for \code{davies} method. 
#' @param tol the eigenvalue cutoff, below which is considered to be 0. This is used to reduce the computation burden.
#' @return
#' \describe{A list containing
#'   \item{p.value}{ association p-value}
#'   \item{Q.adj}{ adjusted score statistic}
#' }
#' @keywords SKAT
#' @importFrom stats as.formula  model.matrix lm resid uniroot pnorm rnorm rbinom binomial
#' @import CompQuadForm
#' @rdname mSKAT2
#' @author Jun Chen
#' @export
#' @references
#' Zhan X, et al. (2017)  A small-sample multivariate kernel machine test for microbiome association studies.
#' Genet Epidemiol. 41(3):210-220.
#' @examples
#' set.seed(123)
#' L <- matrix(rnorm(1000), 100, 10)  # Latent factor to induce correlation in Y
#' Y <- scale(L %*% matrix(rnorm(10 * 40), 10, 40)) + matrix(rnorm(100 * 40), 100, 40)
#' G <- scale(Y %*% matrix(rnorm(40 * 10), 40, 10)) * 0.2 + matrix(rnorm(100 * 10), 100, 10)
#' K <- G %*% t(G)
#' mSKAT2(Y = Y, K = K)  # No de-correlation
#' mSKAT(Y = Y, K = K)   # De-correlation         


mSKAT2 <- function(formula.H0 = NULL, data = NULL, Y, K, acc = 0.00001, lim = 10000, tol = 1e-10) {
	
	n <- ncol(K)
	if (!is.null(formula.H0)) {
		X <- model.matrix(formula.H0, data = data)
		if (nrow(X) == 0) {
			X <- rep(1, n)
		}
	} else {
		X <- rep(1, n)
	}
	
	I <- diag(n)
	P <- I - X %*% solve(t(X) %*% X) %*% t(X)
	Y0 <- P %*% Y
	
	t0 <-  sum(Y0 %*% t(Y0) * K) / sum(Y0 * Y0)
	
	lambda1 <- eigen(P %*% K %*% P - t0 * P, only.values=T)$values
	lambda1 <- lambda1[abs(lambda1) >= tol]
	
	lambda2 <- (svd(Y0 / sqrt(n-1), nu = 0, nv = 0)$d)^2
	lambda2 <- lambda2[lambda2 >= tol]
	lambdas <- as.vector(outer(lambda2, lambda1, '*'))
	p.value <- KAT.pval(0, lambda = sort(lambdas, decreasing=T), acc = acc, lim = lim) 
	
	list(p.value = p.value, Q.adj = t0)
}



#' Small-sample adjusted SKAT for a  correlated (continuous) outcome ('c' stands for 'correlated'). 
#'
#' Compute the adjusted score statistic and p-value for a correlated outcome
#' @param formula.H0  a two-sided linear formula object describing both the fixed-effects and random-effects part of
#'  the model under the null, use the same syntax as the \code{lmer} in \code{lme4} package.
#' @param data  a data frame (required) containing the variables named in formula.
#' @param K    the kernel matrix, which quantifies the similarities between samples.
#' @param relmat   an optional correlation structure, see \code{lme4qtl} package.
#' @param lim maximum number of integration terms for \code{davies} method. 
#' @param acc error bound for \code{davies} method. 
#' @param tol the eigenvalue cutoff, below which it is considered to be 0. This is used to reduce the computation burden.
#' @param ... arguments passing to mixed model \code{glmer} and \code{relmatGlmer}.
#' @return A list containing
#' \describe{ 
#'   \item{p.value}{ association p-value}
#'   \item{Q.adj}{ adjusted score statistic}
#' }
#' @keywords SKAT
#' @importFrom stats as.formula  model.matrix lm resid uniroot pnorm rnorm rbinom binomial
#' @import lme4
#' @import lme4qtl
#' @import CompQuadForm
#' @rdname cSKAT
#' @author Jun Chen
#' @export
#' @references
#' Zhan X, et al. (2018) A small-sample kernel association test for correlated data with application to microbiome association studies.  
#' Genet Epidemiol.;42(8):772-782
#' @examples
#' set.seed(123)
#' Y <- rnorm(100)
#' Z <- matrix(rnorm(200), 100, 2)
#' ID <- gl(20, 5)
#' G <- matrix(rbinom(1000, 1, 0.05), 100, 10)
#' K <- G %*% t(G)
#' data <- data.frame(Y = Y, Z = Z, ID = ID)
#' cSKAT(Y ~ Z + (1 | ID), data = data, K = K)
#' cSKAT(Y ~ (1 | ID), data = data, K = K)  # No covariate


cSKAT <- function (formula.H0, data = NULL, K, relmat = NULL, acc = 0.00001, lim = 10000, tol = 1e-10, ...)  {
	
	H0.lmer <- formula.H0
	temp <- as.character(formula.H0)
	temp[3] <- gsub('\\s*\\+*\\s*\\(.*?\\)\\s*', '',  temp[3], perl = TRUE)
	
	if (temp[3] == '') {
		temp[3] <- '1'
	}
	
	H0.lm <- as.formula(paste(temp[2], temp[1], temp[3]))
	
	if (is.null(relmat)) {
		lmer.obj <- lmer(H0.lmer, data, ...)
	} else {
		lmer.obj <- relmatLmer(H0.lmer, data, relmat = relmat, ...)
	}
	
	var.d <- crossprod(getME(lmer.obj, "Lambdat"))
	Zt <- getME(lmer.obj, "Zt")
	Ztt <- getME(lmer.obj, "Z")
	vr <- sigma(lmer.obj)^2
	
	var.b <- (Ztt %*% var.d %*% Zt)
	sI <- diag(nrow(var.b))
	V2 <- vr * (var.b + sI)
	
	Vi <- sqrt.inv(V2)$Vi
	
	X1 <- model.matrix(H0.lm, data)
	lhs <- H0.lm[[2]]
	Y <- eval(lhs, data)
	
	X1 <- Vi %*% X1
	Y <- Vi %*% Y
	K <- Vi %*% K %*% Vi
	
	mod <- lm(Y ~ X1 - 1)
	res <- resid(mod)
	s2 <- sum(res^2) 
	
	D0  <- diag(length(Y)) 
	P0  <- D0 - X1 %*% solve(t(X1)%*%X1) %*% t(X1)
	PKP <- P0 %*% K %*% P0
	
	q <- as.numeric(res %*% K %*% res / s2)
	ee <- eigen(PKP - q * P0, symmetric = T)  
	lambda0 = ee$values[abs(ee$values) >= tol]
#	p1 <- davies(0, lambda=lambda0)$Qq
#	p.value <- davies(0, lambda=lambda0, acc=1e-9)$Qq
	p.value <- KAT.pval(0, lambda=sort(lambda0, decreasing=T), acc = acc, lim = lim) 
	
	return(list(p.value = p.value, Q.adj = q))
}

#' Small-sample adjusted SKAT for a generalized mixed effects model (GLMM). Since it encompasses a wide range of designs,
#' we name it more generally as Variant-Set Association Test (VSAT).
#' 
#' Compute the adjusted score statistic and p-value for GLMM
#' @param formula.H0  a two-sided linear formula  under the null, using the syntax of the \code{lme4} package.
#' @param data  a data frame (required) containing the variables named in formula.
#' @param K    the kernel matrix, which quantifies the similarities between samples.
#' @param relmat an optional correlation structure, see \code{lme4qtl} package.
#' @param family a GLM family, see \code{glm} and \code{family}. 
#' @param lim maximum number of integration terms for \code{davies} method. 
#' @param acc error bound for \code{davies} method. 
#' @param tol the eigenvalue cutoff, below which it is considered to be 0. This is used to reduce the computation burden.
#' @param ... arguments passing to mixed model \code{glmer} and \code{relmatGlmer}.
#' @return A list containing
#' \describe{
#'   \item{p.value}{association p-value}
#'   \item{Q.adj}{ adjusted score statistic}
#' }
#' @keywords SKAT
#' @importFrom stats as.formula  model.matrix lm resid uniroot pnorm rnorm rbinom binomial predict sigma
#' @importFrom Matrix crossprod 
#' @import lme4
#' @import lme4qtl
#' @import CompQuadForm
#' @rdname VSAT
#' @author Jun Chen
#' @export
#' @references
#' Zhan X, et al. (2020)  Variant-Set Association Test for Generalized Linear Mixed Model.Genet Epidemiol. Submitted.
#' @examples
#' set.seed(123)
#' m <- 20; k <- 4; n <- m * k; s <- 20
#' x <- rnorm(n); z <- matrix(sample(0 : 1, n * s, repl = TRUE, prob = c(0.8, 0.2)), n, s)
#' a <- rnorm(s); b <- 1
#' y <- rbinom(n, 1, prob = binomial()$linkinv(z %*% a + x * b + 0.5 * rnorm(m)[gl(m, k)]))
#' data <- data.frame(y = as.numeric(y), x = x, subject = gl(m, k))
#' rownames(data) <- paste(1:n)
#' K <- z %*% t(z)
#' VSAT(formula.H0 = y ~ x + (1 | subject), data = data, K = K, family = binomial()) 
#' VSAT(formula.H0 = y ~ (1 | subject), data = data, K = K, family = binomial()) # No covariate

VSAT <- function (formula.H0, data = NULL, K, relmat = NULL, family, acc = 0.00001, lim = 10000, tol = 1e-10, ...) {
	
	H0.lmer <- formula.H0
	
	temp <- as.character(formula.H0)
	temp[3] <- gsub('\\s*\\+*\\s*\\(.*?\\)\\s*', '',  temp[3], perl = TRUE)
	
	if (temp[3] == '') {
		temp[3] <- '1'
	}

	H0.lm <- as.formula(paste(temp[2], temp[1], temp[3]))
	
	lhs <- H0.lmer[[2]]
	Y <- eval(lhs, data)
	X1 <- model.matrix(H0.lm, data)
	
	if (is.null(relmat)) {
		lmer.obj <- glmer(H0.lmer, data, family = family, ...)
	} else {
		lmer.obj <- relmatGlmer(H0.lmer, data, relmat = relmat, ...)
	}
	
	eta <- predict(lmer.obj)
	mu <- predict(lmer.obj, type = 'response')
	mu.eta.val <- family$mu.eta(eta)
	Yt <- eta + (Y - mu) / mu.eta.val
	const <- 1
	phi <- 1
	W <- mu.eta.val^2 / family$variance(mu) / const / phi
	
	var.d <- crossprod(getME(lmer.obj, "Lambdat"))
	Zt <- getME(lmer.obj, "Zt")
	Ztt <- getME(lmer.obj, "Z")
	
	var.b <- (Ztt %*% var.d %*% Zt)
	V2 <- (var.b + diag(1 / W))
	
	eig.obj <- eigen(V2, symmetric = TRUE)
	vectors <- eig.obj$vectors
	values <- eig.obj$values
	ind <- values >= tol
	values <- values[ind]
	vectors <- vectors[, ind]
	temp <- t(vectors) * sqrt(values)
	V <- vectors  %*% temp
	
	temp <- t(vectors) / sqrt(values)
	Vi2 <- t(temp)  %*% temp
	Vi <- vectors  %*% temp
	
	P <- Vi2 - Vi2 %*% X1 %*% solve(t(X1) %*% Vi2 %*% X1) %*% t(X1)  %*% Vi2
	PKP <- P %*% K %*% P 
	PP <-  P  %*% P
	
	VPKPV <- V %*% PKP %*% V
	VPPV <- V %*% PP %*% V
	
	q <- as.vector((t(Yt) %*% PKP %*% Yt) / (t(Yt) %*% PP %*% Yt))
	
	ee <- eigen(VPKPV - q * VPPV, symmetric = T) 
	
	lambda <- ee$values[abs(ee$values) >= tol]
	#p.value <- davies(0, lambda=lambda)$Qq
	p.value <- KAT.pval(0, lambda = sort(lambda, decreasing=T), acc = acc, lim = lim) 
	return(list(p.value = p.value,  Q.adj = q))
	
}

