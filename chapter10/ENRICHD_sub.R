# Functions for the Sub models
# Authored by: Dr. Wenhua Jiang
# Reference: Wu, C. O., Tian, X. and Jiang, W. A shared parameter model for the
#   estimation of longitudinal concomitant intervention effects Biostatistics,
#   12(4):737–749, 2011.
#######################################################
## apply function to each row of  Samples
#######################################################

simcomp <- function(like.fun, samples, ...) {
	likelihood <- apply(samples, 1, like.fun, ...)
	return(likelihood)
}

#####################################################
## conditional likelihood of y given a, b, sigma.y ##
#####################################################

y.likelihood.ab <- function(ab, y, Z.0, Z.1, sigma.y, lambda) {
	n <- length(y)
	a <- ab[1:dim(Z.0)[2]]
	b <- ab[(dim(Z.0)[2] + 1):(dim(Z.0)[2] + dim(Z.1)[2])]
	resid <- y - Z.0 %*% a - lambda * (Z.1 %*% b)
	fy.ab <- (1/(((2 * pi)^(n/2)) * ((sigma.y^2)^(n/2)))) * exp(-(1/(2 * sigma.y^2)) * sum(resid^2))
	return(fy.ab)
}

y.likelihood.a <- function(a, y, Z.0, sigma.y) {
	n <- length(y)
	resid <- y - Z.0 %*% a
	fy.a <- (1/(((2 * pi)^(n/2)) * ((sigma.y^2)^(n/2)))) * exp(-(1/(2 * sigma.y^2)) * sum(resid^2))
	return(fy.a)
}

#################################################################
## compute conditional likelihood of S given gamma, a, sigma.s ##
#################################################################

S.likelihood.0 <- function(U, gamma, sigma.s, S) {
	resid <- S - gamma %*% U
	fS.0 <- (1/(((2 * pi)^(1/2)) * ((sigma.s^2)^(1/2)))) * exp(-(1/(2 * sigma.s^2)) * (resid^2))
	return(fS.0)
}

S.likelihood.1 <- function(U, gamma, sigma.s, t.point) {
	t.max <- max(t.point)
	mu <- gamma %*% U
	fS.1 <- 1 - pnorm(t.max, mu, sigma.s)
	return(fS.1)
}

S.likelihood.2 <- function(U, gamma, sigma.s, t.point) {
	t.min <- min(t.point)
	mu <- gamma %*% U
	fS.2 <- pnorm(t.min, mu, sigma.s)
	return(fS.2)
}

#########################################
## inner integration with respect to S ##
#########################################

g.int.S.1 <- function(U, gamma, sigma.s, t.point) {
	t.max <- max(t.point)
	g.1 <- (1/(((2 * pi)^(1/2)) * ((sigma.s^2)^(1/2)))) * exp(-(t.max - gamma %*% U)^2/(2 * sigma.s^2))
	return(g.1)
}

g.int.S.2 <- function(U, gamma, sigma.s, t.point) {
	t.max <- max(t.point)
	g.2 <- (1/(((2 * pi)^(1/2)) * ((sigma.s^2)^(1/2)))) * exp(-(t.max - gamma %*% U)^2/(2 * sigma.s^2)) * ((t.max - gamma %*% U)/sigma.s^2)
	return(g.2)
}

g.int.S.3 <- function(U, gamma, sigma.s, t.point) {
	t.min <- min(t.point)
	g.3 <- -(1/(((2 * pi)^(1/2)) * ((sigma.s^2)^(1/2)))) * exp(-(t.min - gamma %*% U)^2/(2 * sigma.s^2))
	return(g.3)
}

g.int.S.4 <- function(U, gamma, sigma.s, t.point) {
	t.min <- min(t.point)
	g.4 <- -(1/(((2 * pi)^(1/2)) * ((sigma.s^2)^(1/2)))) * exp(-(t.min - gamma %*% U)^2/(2 * sigma.s^2)) * ((t.min - gamma %*% U)/sigma.s^2)
	return(g.4)
}

######################################################
## sample a and b, compute all the necessary vector ##
######################################################

generate.ab <- function(n.sample, alpha, beta, cov.a, cov.b) {
    ## revise here when change model
	beta <- c(beta[1], 0, 0, 0, 0, beta[2])
	beta <- matrix(c(beta), nrow = 2, ncol = 3, byrow = T)
	mean.a <- alpha
	mean.b <- c(beta %*% c(1, alpha))
	mean.ab <- c(mean.a, mean.b)
	cov.ab <- matrix(NA, dim(cov.a)[1] + dim(cov.b)[1], dim(cov.a)[2] + dim(cov.b)[2])
	cov.ab[1:dim(cov.a)[1], 1:dim(cov.a)[2]] <- cov.a
	cov.ab[(dim(cov.a)[1] + 1):dim(cov.ab)[1], (dim(cov.a)[2] + 1):dim(cov.ab)[2]] <- cov.b + beta %*% rbind(0, cbind(0, cov.a)) %*% t(beta)
	cov.ab[1:dim(cov.a)[1], (dim(cov.a)[2] + 1):dim(cov.ab)[2]] <- cbind(0, cov.a) %*% t(beta)
	cov.ab[(dim(cov.a)[1] + 1):dim(cov.ab)[1], 1:dim(cov.a)[2]] <- beta %*% rbind(0, cov.a)
	ab <- mvrnorm(n.sample, mu = mean.ab, Sigma = cov.ab)
	a <- ab[, 1:length(mean.a)]
	b <- ab[, (length(mean.a) + 1):length(mean.ab)]
	U <- cbind(1, a[, 1])
	sample.ab = list(a = a, b = b, U = U)
	return(sample.ab)
}

##################################
## compute useful quantitatives ##
##################################

generate.quant <- function(data, sample.ab, alpha, beta, gamma, sigma.s, sigma.y) {
	a <- sample.ab$a
	b <- sample.ab$b
	ab <- cbind(a, b)
	U <- sample.ab$U
	S <- data[1, 8]
	n.mes <- data[2, 8]
	t.point <- data[1:n.mes, 1]
	y <- data[1:n.mes, 2]
	Z.0 <- data[1:n.mes, 3:4]
	Z.1 <- data[1:n.mes, 5:6]
	lambda <- data[1:n.mes, 7]
	fy.ab <- simcomp("y.likelihood.ab", samples = ab, y = y, Z.0 = Z.0, Z.1 = Z.1, sigma.y = sigma.y, lambda = lambda)
	fy.a <- simcomp("y.likelihood.a", samples = a, y = y, Z.0 = Z.0, sigma.y = sigma.y)
	fS.0 <- simcomp("S.likelihood.0", samples = U, gamma = gamma, sigma.s = sigma.s, S = S)
	fS.1 <- simcomp("S.likelihood.1", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	fS.2 <- simcomp("S.likelihood.2", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	g.1 <- simcomp("g.int.S.1", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	g.2 <- simcomp("g.int.S.2", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	g.3 <- simcomp("g.int.S.3", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	g.4 <- simcomp("g.int.S.4", samples = U, gamma = gamma, sigma.s = sigma.s, t.point = t.point)
	factor.gamma <- (as.vector(S - U %*% gamma) * U)/sigma.s^2
	data.obj = list(fy.ab = fy.ab, fy.a = fy.a, fS.0 = fS.0, fS.1 = fS.1, fS.2 = fS.2, g.1 = g.1, g.2 = g.2, g.3 = g.3, g.4 = g.4, factor.gamma = factor.gamma)
	return(data.obj)
}

generate.factor <- function(sample.ab, alpha, beta, cov.a, cov.b) {
    ## revise here when change model
	beta <- c(beta[1], 0, 0, 0, 0, beta[2])
	beta <- matrix(c(beta), nrow = 2, ncol = 3, byrow = T)
	a <- sample.ab$a
	b <- sample.ab$b
	factor.alpha <- t(solve(cov.a) %*% (t(a) - alpha))
    ## revise here when change model
	V.a.1 <- rep(1, dim(a)[1])
	V.a.2 <- a[, 2]
	V.a <- cbind(1, a)
	factor.beta.coeff <- t(solve(cov.b) %*% (t(b) - beta %*% t(V.a)))
	factor.beta <- cbind(factor.beta.coeff[, 1] * V.a.1, factor.beta.coeff[, 2] * V.a.2)
	data.factor = list(factor.alpha = factor.alpha, factor.beta = factor.beta)
	return(data.factor)
}

###################################
## compute the type I likelihood ##
###################################

likelihood.type.1 <- function(data.obj) {
	fy.ab <- data.obj$fy.ab
	fS.0 <- data.obj$fS.0
	likelihood <- mean(fy.ab * fS.0)
	return(likelihood)
}

###################################################
## compute the gradient of the type I likelihood ##
###################################################

gradient.type.1 <- function(data.obj, data.factor, sample.ab) {
	a <- sample.ab$a
	b <- sample.ab$b
	fy.ab <- data.obj$fy.ab
	fS.0 <- data.obj$fS.0
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	likelihood <- mean(fy.ab * fS.0)
	partial.alpha <- apply(fy.ab * fS.0 * factor.alpha, 2, mean)
	gradient.alpha <- partial.alpha/likelihood
	partial.beta <- apply(fy.ab * fS.0 * factor.beta, 2, mean)
	gradient.beta <- partial.beta/likelihood
	partial.gamma <- apply(fy.ab * fS.0 * factor.gamma, 2, mean)
	gradient.gamma <- as.vector(partial.gamma/likelihood)
	gradient <- c(gradient.alpha, gradient.beta, gradient.gamma)
	return(gradient)
}

##################################################
## compute the Hessian of the type I likelihood ##
##################################################

Hessian.type.1 <- function(data.obj, data.factor, gradient, sample.ab, cov.a, cov.b, sigma.s) {
	a <- sample.ab$a
	b <- sample.ab$b
	U <- sample.ab$U
	fy.ab <- data.obj$fy.ab
	fS.0 <- data.obj$fS.0
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	l.alpha <- dim(factor.alpha)[2]
	l.beta <- dim(factor.beta)[2]
	l.gamma <- dim(factor.gamma)[2]
	gradient.alpha <- gradient[1:l.alpha]
	gradient.beta <- gradient[(l.alpha + 1):(l.alpha + l.beta)]
	gradient.gamma <- gradient[(l.alpha + l.beta + 1):length(gradient)]
	n.sample <- dim(a)[1]
	Hessian <- matrix(NA, length(gradient), length(gradient))
	likelihood <- mean(fy.ab * fS.0)	
	factor <- fy.ab * fS.0
	factor.alpha.alpha <- factor.alpha * sqrt(factor)
	partial.alpha.alpha <- (t(factor.alpha.alpha) %*% factor.alpha.alpha)/n.sample - mean(factor) * solve(cov.a)
	Hessian.alpha <- partial.alpha.alpha/likelihood - (gradient.alpha %o% gradient.alpha)
	factor.beta.beta.1 <- factor.beta * sqrt(factor)
	inv.cov.b <- solve(cov.b)
    ## revise here when change model
	factor.beta.beta.2 <- matrix(NA, 2, 2)
	factor.beta.beta.2[1, 1] <- inv.cov.b[1, 1]
	factor.beta.beta.2[2, 2] <- inv.cov.b[2, 2]
	factor.beta.beta.2[1, 2] <- inv.cov.b[1, 2]
	factor.beta.beta.2[2, 1] <- inv.cov.b[2, 1]
	V.a.1 <- rep(1, dim(a)[1])
	V.a.2 <- a[, 2]
	V.a <- cbind(V.a.1, V.a.2)
	factor.beta.beta.3 <- sqrt(factor) * V.a
	partial.beta.beta <- (t(factor.beta.beta.1) %*% factor.beta.beta.1)/n.sample - (t(factor.beta.beta.3) %*% factor.beta.beta.3)/n.sample * factor.beta.beta.2
	Hessian.beta <- partial.beta.beta/likelihood - (gradient.beta %o% gradient.beta)
	factor.gamma.gamma.1 <- factor.gamma * sqrt(factor)
	factor.gamma.gamma.2 <- (sqrt(factor) * U)/sigma.s
	partial.gamma.gamma <- (t(factor.gamma.gamma.1) %*% factor.gamma.gamma.1)/n.sample - (t(factor.gamma.gamma.2) %*% factor.gamma.gamma.2)/n.sample
	Hessian.gamma <- partial.gamma.gamma/likelihood - (gradient.gamma %o% gradient.gamma)
	factor.alpha.beta.1 <- factor.alpha * sqrt(factor)
	factor.alpha.beta.2 <- factor.beta * sqrt(factor)
	partial.alpha.beta <- (t(factor.alpha.beta.1) %*% factor.alpha.beta.2)/n.sample
	Hessian.alpha.beta <- partial.alpha.beta/likelihood - (gradient.alpha %o% gradient.beta)
	factor.alpha.gamma.1 <- factor.alpha * sqrt(factor)
	factor.alpha.gamma.2 <- factor.gamma * sqrt(factor)
	partial.alpha.gamma <- (t(factor.alpha.gamma.1) %*% factor.alpha.gamma.2)/n.sample
	Hessian.alpha.gamma <- partial.alpha.gamma/likelihood - (gradient.alpha %o% gradient.gamma)
	factor.beta.gamma.1 <- factor.beta * sqrt(factor)
	factor.beta.gamma.2 <- factor.gamma * sqrt(factor)
	partial.beta.gamma <- (t(factor.beta.gamma.1) %*% factor.beta.gamma.2)/n.sample
	Hessian.beta.gamma <- partial.beta.gamma/likelihood - (gradient.beta %o% gradient.gamma)
	Hessian[1:l.alpha, 1:l.alpha] <- Hessian.alpha
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.beta
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.gamma
	Hessian[1:l.alpha, (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.alpha.beta
	Hessian[1:l.alpha, (l.alpha + l.beta + 1):length(gradient)] <- Hessian.alpha.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.beta.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), 1:l.alpha] <- t(Hessian.alpha.beta)
	Hessian[(l.alpha + l.beta + 1):length(gradient), 1:l.alpha] <- t(Hessian.alpha.gamma)
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + 1):(l.alpha + l.beta)] <- t(Hessian.beta.gamma)
	return(Hessian)
}

####################################
## compute the type II likelihood ##
####################################

likelihood.type.2 <- function(data.obj) {
	fy.a <- data.obj$fy.a
	fS.1 <- data.obj$fS.1
	likelihood <- mean(fy.a * fS.1)
	return(likelihood)
}

####################################################
## compute the gradient of the type II likelihood ##
####################################################

gradient.type.2 <- function(data.obj, data.factor, sample.ab) {
	a <- sample.ab$a
	U <- sample.ab$U
	fy.a <- data.obj$fy.a
	fS.1 <- data.obj$fS.1
	g.1 <- data.obj$g.1
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	likelihood <- mean(fy.a * fS.1)
	partial.alpha <- apply(fy.a * fS.1 * factor.alpha, 2, mean)
	gradient.alpha <- partial.alpha/likelihood
	gradient.beta <- rep(0, dim(factor.beta)[2])
	partial.gamma <- apply(fy.a * g.1 * U, 2, mean)
	gradient.gamma <- as.vector(partial.gamma/likelihood)
	gradient <- c(gradient.alpha, gradient.beta, gradient.gamma)
	return(gradient)
}

###################################################
## compute the Hessian of the type II likelihood ##
###################################################

Hessian.type.2 <- function(data.obj, data.factor, gradient, sample.ab, cov.a, cov.b, sigma.s) {
	a <- sample.ab$a
	U <- sample.ab$U
	fy.a <- data.obj$fy.a
	fS.1 <- data.obj$fS.1
	g.1 <- data.obj$g.1
	g.2 <- data.obj$g.2
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	l.alpha <- dim(factor.alpha)[2]
	l.beta <- dim(factor.beta)[2]
	l.gamma <- dim(factor.gamma)[2]
	gradient.alpha <- gradient[1:l.alpha]
	gradient.beta <- gradient[(l.alpha + 1):(l.alpha + l.beta)]
	gradient.gamma <- gradient[(l.alpha + l.beta + 1):length(gradient)]
	n.sample <- dim(a)[1]
	Hessian <- matrix(NA, length(gradient), length(gradient))
	likelihood <- mean(fy.a * fS.1)
	factor <- fy.a * fS.1
	factor.alpha.alpha <- factor.alpha * sqrt(factor)
	partial.alpha.alpha <- (t(factor.alpha.alpha) %*% factor.alpha.alpha)/n.sample - mean(factor) * solve(cov.a)
	Hessian.alpha <- partial.alpha.alpha/likelihood - (gradient.alpha %o% gradient.alpha)
	Hessian.beta <- matrix(0, l.beta, l.beta)
	factor.gamma.gamma.1 <- (fy.a * g.2) * U
	factor.gamma.gamma.2 <- U
	partial.gamma.gamma <- (t(factor.gamma.gamma.1) %*% factor.gamma.gamma.2)/n.sample
	Hessian.gamma <- partial.gamma.gamma/likelihood - (gradient.gamma %o% gradient.gamma)
	Hessian.alpha.beta <- matrix(0, l.alpha, l.beta)
	factor.alpha.gamma.1 <- factor.alpha
	factor.alpha.gamma.2 <- (fy.a * g.1) * U
	partial.alpha.gamma <- (t(factor.alpha.gamma.1) %*% factor.alpha.gamma.2)/n.sample
	Hessian.alpha.gamma <- partial.alpha.gamma/likelihood - (gradient.alpha %o% gradient.gamma)
	Hessian.beta.gamma <- matrix(0, l.beta, l.gamma)
	Hessian[1:l.alpha, 1:l.alpha] <- Hessian.alpha
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.beta
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.gamma
	Hessian[1:l.alpha, (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.alpha.beta
	Hessian[1:l.alpha, (l.alpha + l.beta + 1):length(gradient)] <- Hessian.alpha.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.beta.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), 1:l.alpha] <- t(Hessian.alpha.beta)
	Hessian[(l.alpha + l.beta + 1):length(gradient), 1:l.alpha] <- t(Hessian.alpha.gamma)
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + 1):(l.alpha + l.beta)] <- t(Hessian.beta.gamma)
	return(Hessian)
}

#####################################
## compute the type III likelihood ##
#####################################

likelihood.type.3 <- function(data.obj) {
	fy.ab <- data.obj$fy.ab
	fS.2 <- data.obj$fS.2
	likelihood <- mean(fy.ab * fS.2)
	return(likelihood)
}

#####################################################
## compute the gradient of the type III likelihood ##
#####################################################

gradient.type.3 <- function(data.obj, data.factor, sample.ab) {
	a <- sample.ab$a
	b <- sample.ab$b
	U <- sample.ab$U
	fy.ab <- data.obj$fy.ab
	fS.2 <- data.obj$fS.2
	g.3 <- data.obj$g.3
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	likelihood <- mean(fy.ab * fS.2)
	partial.alpha <- apply(fy.ab * fS.2 * factor.alpha, 2, mean)
	gradient.alpha <- partial.alpha/likelihood	
	partial.beta <- apply(fy.ab * fS.2 * factor.beta, 2, mean)
	gradient.beta <- partial.beta/likelihood
	partial.gamma <- apply(fy.ab * g.3 * U, 2, mean)
	gradient.gamma <- as.vector(partial.gamma/likelihood)
	gradient <- c(gradient.alpha, gradient.beta, gradient.gamma)
	return(gradient)
}

####################################################
## compute the Hessian of the type III likelihood ##
####################################################

Hessian.type.3 <- function(data.obj, data.factor, gradient, sample.ab, cov.a, cov.b, sigma.s) {
	a <- sample.ab$a
	b <- sample.ab$b
	U <- sample.ab$U
	fy.ab <- data.obj$fy.ab
	fS.2 <- data.obj$fS.2
	g.3 <- data.obj$g.3
	g.4 <- data.obj$g.4
	factor.alpha <- data.factor$factor.alpha
	factor.beta <- data.factor$factor.beta
	factor.gamma <- data.obj$factor.gamma
	l.alpha <- dim(factor.alpha)[2]
	l.beta <- dim(factor.beta)[2]
	l.gamma <- dim(factor.gamma)[2]
	gradient.alpha <- gradient[1:l.alpha]
	gradient.beta <- gradient[(l.alpha + 1):(l.alpha + l.beta)]
	gradient.gamma <- gradient[(l.alpha + l.beta + 1):length(gradient)]
	n.sample <- dim(a)[1]
	Hessian <- matrix(NA, length(gradient), length(gradient))
	likelihood <- mean(fy.ab * fS.2)
	factor <- fy.ab * fS.2
	factor.alpha.alpha <- factor.alpha * sqrt(factor)
	partial.alpha.alpha <- (t(factor.alpha.alpha) %*% factor.alpha.alpha)/n.sample - mean(factor) * solve(cov.a)
	Hessian.alpha <- partial.alpha.alpha/likelihood - (gradient.alpha %o% gradient.alpha)
	factor.beta.beta.1 <- factor.beta * sqrt(factor)
	inv.cov.b <- solve(cov.b)
    ## revise here when change model
	factor.beta.beta.2 <- matrix(NA, 2, 2)
	factor.beta.beta.2[1, 1] <- inv.cov.b[1, 1]
	factor.beta.beta.2[2, 2] <- inv.cov.b[2, 2]
	factor.beta.beta.2[1, 2] <- inv.cov.b[1, 2]
	factor.beta.beta.2[2, 1] <- inv.cov.b[2, 1]
	V.a.1 <- rep(1, dim(a)[1])
	V.a.2 <- a[, 2]
	V.a <- cbind(V.a.1, V.a.2)
	factor.beta.beta.3 <- sqrt(factor) * V.a
	partial.beta.beta <- (t(factor.beta.beta.1) %*% factor.beta.beta.1)/n.sample - (t(factor.beta.beta.3) %*% factor.beta.beta.3)/n.sample * factor.beta.beta.2
	Hessian.beta <- partial.beta.beta/likelihood - (gradient.beta %o% gradient.beta)
	factor.gamma.gamma.1 <- (fy.ab * g.4) * U
	factor.gamma.gamma.2 <- U
	partial.gamma.gamma <- (t(factor.gamma.gamma.1) %*% factor.gamma.gamma.2)/n.sample
	Hessian.gamma <- partial.gamma.gamma/likelihood - (gradient.gamma %o% gradient.gamma)
	factor.alpha.beta.1 <- factor.alpha * sqrt(factor)
	factor.alpha.beta.2 <- factor.beta * sqrt(factor)
	partial.alpha.beta <- (t(factor.alpha.beta.1) %*% factor.alpha.beta.2)/n.sample
	Hessian.alpha.beta <- partial.alpha.beta/likelihood - (gradient.alpha %o% gradient.beta)
	factor.alpha.gamma.1 <- factor.alpha
	factor.alpha.gamma.2 <- (fy.ab * g.3) * U
	partial.alpha.gamma <- (t(factor.alpha.gamma.1) %*% factor.alpha.gamma.2)/n.sample
	Hessian.alpha.gamma <- partial.alpha.gamma/likelihood - (gradient.alpha %o% gradient.gamma)
	factor.beta.gamma.1 <- factor.beta
	factor.beta.gamma.2 <- (fy.ab * g.3) * U
	partial.beta.gamma <- (t(factor.beta.gamma.1) %*% factor.beta.gamma.2)/n.sample
	Hessian.beta.gamma <- partial.beta.gamma/likelihood - (gradient.beta %o% gradient.gamma)
	Hessian[1:l.alpha, 1:l.alpha] <- Hessian.alpha
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.beta
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.gamma
	Hessian[1:l.alpha, (l.alpha + 1):(l.alpha + l.beta)] <- Hessian.alpha.beta
	Hessian[1:l.alpha, (l.alpha + l.beta + 1):length(gradient)] <- Hessian.alpha.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), (l.alpha + l.beta + 1):length(gradient)] <- Hessian.beta.gamma
	Hessian[(l.alpha + 1):(l.alpha + l.beta), 1:l.alpha] <- t(Hessian.alpha.beta)
	Hessian[(l.alpha + l.beta + 1):length(gradient), 1:l.alpha] <- t(Hessian.alpha.gamma)
	Hessian[(l.alpha + l.beta + 1):length(gradient), (l.alpha + 1):(l.alpha + l.beta)] <- t(Hessian.beta.gamma)
	return(Hessian)
}

###################################################
## compute the log-likelihood for single subject ##
###################################################

log.likelihood.sub <- function(data, sample.ab, alpha, beta, gamma, sigma.s, sigma.y) {
	S <- data[1, 8]
	n.mes <- data[2, 8]
	t.point <- data[1:n.mes, 1]
	data.obj <- generate.quant(data = data, sample.ab = sample.ab, alpha = alpha, beta = beta, gamma = gamma, sigma.s = sigma.s, sigma.y = sigma.y)
	if (S >= min(t.point) && S <= max(t.point)) {
		log.likelihood <- log(likelihood.type.1(data.obj = data.obj))
	}
	if (S > max(t.point)) {
		log.likelihood <- log(likelihood.type.2(data.obj = data.obj))
	}
	if (S < min(t.point)) {
		log.likelihood <- log(likelihood.type.3(data.obj = data.obj))
	}
	result <- log.likelihood
	return(result)
}

################################
## compute the log-likelihood ##
################################

log.likelihood <- function(data, n.sample, alpha, beta, gamma, cov.a, cov.b, sigma.s, sigma.y) {
	sample.ab <- generate.ab(n.sample = n.sample, alpha = alpha, beta = beta, cov.a = cov.a, cov.b = cov.b)
	result <- apply(data, 3, log.likelihood.sub, sample.ab = sample.ab, alpha = alpha, beta = beta, gamma = gamma, sigma.s = sigma.s, sigma.y = sigma.y)
	log.likelihood <- sum(result)
	return(log.likelihood)
}

######################################################
## compute the gradient, Hessian for single subject ##
######################################################

gradient.Hessian.sub <- function(data, data.factor, sample.ab, alpha, beta, gamma, cov.a, cov.b, sigma.s, sigma.y) {
	S <- data[1, 8]
	n.mes <- data[2, 8]
	t.point <- data[1:n.mes, 1]
	data.obj <- generate.quant(data = data, sample.ab = sample.ab, alpha = alpha, beta = beta, gamma = gamma, sigma.s = sigma.s, sigma.y = sigma.y)
	if (S >= min(t.point) && S <= max(t.point)) {
		log.likelihood <- log(likelihood.type.1(data.obj = data.obj))
		gradient <- gradient.type.1(data.obj = data.obj, data.factor = data.factor, sample.ab = sample.ab)
		Hessian <- Hessian.type.1(data.obj = data.obj, data.factor = data.factor, gradient = gradient, sample.ab = sample.ab, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s)
	}
	if (S > max(t.point)) {
		log.likelihood <- log(likelihood.type.2(data.obj = data.obj))
		gradient <- gradient.type.2(data.obj = data.obj, data.factor = data.factor, sample.ab = sample.ab)
		Hessian <- Hessian.type.2(data.obj = data.obj, data.factor = data.factor, gradient = gradient, sample.ab = sample.ab, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s)
	}
	if (S < min(t.point)) {
		log.likelihood <- log(likelihood.type.3(data.obj = data.obj))
		gradient <- gradient.type.3(data.obj = data.obj, data.factor = data.factor, sample.ab = sample.ab)
		Hessian <- Hessian.type.3(data.obj = data.obj, data.factor = data.factor, gradient = gradient, sample.ab = sample.ab, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s)
	}
	result <- cbind(gradient, Hessian, c(log.likelihood, rep(NA, length(gradient) - 1)))
	dimnames(result) <- list(NULL, NULL)
	return(result)
}

###################################
## compute the gradient, Hessian ##
###################################

gradient.Hessian <- function(data, n.sample, alpha, beta, gamma, cov.a, cov.b, sigma.s, sigma.y) {
	n.par <- length(alpha) + length(beta) + length(gamma)
	sample.ab <- generate.ab(n.sample = n.sample, alpha = alpha, beta = beta, cov.a = cov.a, cov.b = cov.b)
	data.factor <- generate.factor(sample.ab = sample.ab, alpha = alpha, beta = beta, cov.a = cov.a, cov.b = cov.b)
	result <- apply(data, 3, gradient.Hessian.sub, data.factor = data.factor, sample.ab = sample.ab, alpha = alpha, beta = beta, gamma = gamma, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s, sigma.y = sigma.y)
	result <- apply(result, 1, sum)
	gradient <- c(result[1:n.par])
	Hessian <- matrix(result[(n.par + 1):(n.par * n.par + n.par)], n.par, n.par)
	log.likelihood <- result[n.par * n.par + n.par + 1]
	obj <- list(gradient = gradient, Hessian = Hessian, log.likelihood = log.likelihood)
	return(obj)
}

##############################################
## Newton-Raphson algorithm to seek the MLE ##
##############################################

NR.share.parameter <- function(data, n.iter, n.sample, alpha.int, beta.int, gamma.int, cov.a, cov.b, sigma.s, sigma.y, max.loglik = T, trace = T) {
	result <- matrix(NA, n.iter + 1, length(alpha.int) + length(beta.int) + length(gamma.int) + 1)
	alpha.item <- paste(rep("alpha", length(alpha.int)), as.character(0:(length(alpha.int) - 1)))
    ## revise here when change model
	beta.item <- paste(rep("beta", length(beta.int)), c("00", "12"))
	gamma.item <- paste(rep("gamma", length(gamma.int)), as.character(0:(length(gamma.int) - 1)))
	item <- c("log-likelihood", alpha.item, beta.item, gamma.item)
	dimnames(result) <- list(NULL, item)
	alpha.cur <- alpha.int
	beta.cur <- beta.int
	gamma.cur <- gamma.int
	normal <- T
	for (i in 1:n.iter) {
		result[i, 2:(length(alpha.cur) + 1)] <- c(alpha.cur)
		result[i, (length(alpha.cur) + 2):(length(alpha.cur) + length(beta.cur) + 1)] <- c(beta.cur)
		result[i, (length(alpha.cur) + length(beta.cur) + 2):dim(result)[2]] <- c(gamma.cur)
		obj <- gradient.Hessian(data = data, n.sample = n.sample, alpha = alpha.cur, beta = beta.cur, gamma = gamma.cur, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s, sigma.y = sigma.y)
		result[i, 1] <- obj$log.likelihood
		if (trace) {
			cat(i - 1, "log.likelihood", round(result[i, 1], 3), "alpha", round(alpha.cur, 3), "beta", round(beta.cur, 3), "gamma", round(gamma.cur, 3), "\n")
		}
		if (result[i, 1] == -Inf) {
			normal <- F
			break
		}
		eigen.obj <- eigen(obj$Hessian)
		change <- c((eigen.obj$vectors) %*% diag(c(1/eigen.obj$values)) %*% t(eigen.obj$vectors) %*% c(obj$gradient))
		alpha.cur <- alpha.cur - change[1:length(alpha.cur)]
		beta.cur <- beta.cur - change[(length(alpha.cur) + 1):(length(alpha.cur) + length(beta.cur))]
		gamma.cur <- gamma.cur - change[(length(alpha.cur) + length(beta.cur) + 1):length(change)]
		if (max(abs(change)) < 1e-3) break
	}
	if (normal) {
		result[i + 1, 2:(length(alpha.cur) + 1)] <- c(alpha.cur)
		result[i + 1, (length(alpha.cur) + 2):(length(alpha.cur) + length(beta.cur) + 1)] <- c(beta.cur)
		result[i + 1, (length(alpha.cur) + length(beta.cur) + 2):dim(result)[2]] <- c(gamma.cur)
		if (max.loglik) {
			result[i + 1, 1] <- log.likelihood(data = data, n.sample = n.sample, alpha = alpha.cur, beta = beta.cur, gamma = gamma.cur, cov.a = cov.a, cov.b = cov.b, sigma.s = sigma.s, sigma.y = sigma.y)
			if (trace) {
				cat(i, "log.likelihood", round(result[i + 1, 1], 3), "alpha", round(alpha.cur, 3), "beta", round(beta.cur, 3), "gamma", round(gamma.cur, 3), "\n")
			}
			result <- result[1:(i + 1), ]
			parameter <- c(result[c(result[, 1]) == max(c(result[, 1])), ])
		}
		else {
			if (trace) {
				cat(i, "log.likelihood", round(result[i + 1, 1], 3), "alpha", round(alpha.cur, 3), "beta", round(beta.cur, 3), "gamma", round(gamma.cur, 3), "\n")
			}
			result <- result[1:(i + 1), ]
			parameter <- c(result[dim(result)[1], ])
		}
	}
	if (!normal) {
		if (max.loglik) {
			result <- result[1:i, ]
			parameter <- c(result[c(result[, 1]) == max(c(result[, 1])), ])
		}
		else {
			result <- result[1:i, ]
			parameter <- c(result[dim(result)[1], ])
		}
	}
	log.likelihood <- parameter[1]
	alpha <- parameter[2:(length(alpha.cur) + 1)]
	beta <- parameter[(length(alpha.cur) + 2):(length(alpha.cur) + length(beta.cur) + 1)]
	gamma <- parameter[(length(alpha.cur) + length(beta.cur) + 2):dim(result)[2]]
	if (trace) {
		cat("*", "log.likelihood", round(log.likelihood, 3), "alpha", round(alpha, 3), "beta", round(beta, 3), "gamma", round(gamma, 3), "\n")
	}
	whole <- result
	ultimate <- parameter
	obj <- list(whole = whole, ultimate = ultimate)
	return(obj)
}

####################################################
## convert the data form to that required by REML ##
####################################################

convert.data <- function(data) {
	n.mes <- dim(data)[1]
	n.sub <- dim(data)[3]
	t.point <- rep(NA, n.mes * n.sub)
	y <- rep(NA, n.mes * n.sub)
	R <- rep(NA, n.mes * n.sub)
	S <- rep(NA, n.mes * n.sub)
	lambda <- rep(NA, n.mes * n.sub)
	subject <- rep(NA, n.mes * n.sub)
	ind.start <- 0
	ind.end <- 0
	for (i in 1:n.sub) {
		n.obs <- data[2, 8, i]
		ind.start <- ind.end + 1
		ind.end <- ind.end + n.obs
		t.point[ind.start:ind.end] <- data[1:n.obs, 1, i]
		y[ind.start:ind.end] <- data[1:n.obs, 2, i]
		R[ind.start:ind.end] <- data[1:n.obs, 6, i]
		S[ind.start:ind.end] <- rep(data[1, 8, i], n.obs)
		lambda[ind.start:ind.end] <- data[1:n.obs, 7, i]
		subject[ind.start:ind.end] <- rep(i, n.obs)
	}
	lambda.R <- lambda * R
	t.point <- t.point[1:ind.end]
	y <- y[1:ind.end]
	R <- R[1:ind.end]
	S <- S[1:ind.end]
	lambda <- lambda[1:ind.end]
	subject <- subject[1:ind.end]
	lambda.R <- lambda.R[1:ind.end]
	mydat.lme <- groupedData(y ~ t.point + lambda + lambda.R | subject, data = data.frame(y, t.point, lambda, lambda.R, R, S, subject))
	return(mydat.lme)
}

####################################################################################
## fit shared-parameter model by using the initial from linear mixed effect model ##
####################################################################################

fit.lme.sp <- function(data, n.iter, n.sample, max.loglik = T, trace = T) {
	initial <- F
	start.min <- min(data[2, 8, ])
	while (initial == F) {
		data.sp <- data[, , data[2, 8, ] >= start.min]
		data.lme <- convert.data(data.sp)
		fm.lme <- try(lme(fixed = y ~ t.point + lambda + lambda.R, data = data.lme, method = "REML"), silent = T)
		if (inherits(fm.lme, "try-error") == F) initial <- T
		start.min <- start.min + 1
	}
	alpha.lme.int <- fm.lme$coeff$fixed[1:2]
	beta.lme.int <- fm.lme$coeff$fixed[3:4]
	alpha.int <- alpha.lme.int
	cov.hat <- VarCorr(fm.lme)
	cov.a.hat <- diag(c(as.numeric(cov.hat[1:2, 1])))
	cov.a.hat[2, 1] <- sqrt(cov.a.hat[1, 1]) * sqrt(cov.a.hat[2, 2]) * as.numeric(cov.hat[2, 3])
	cov.a.hat[1, 2] <- cov.a.hat[2, 1]
	cov.b.hat <- diag(c(as.numeric(cov.hat[3:4, 1])))
	cov.b.hat[2, 1] <- sqrt(cov.b.hat[1, 1]) * sqrt(cov.b.hat[2, 2]) * as.numeric(cov.hat[4, 5])
	cov.b.hat[1, 2] <- cov.b.hat[2, 1]
	sigma.y.hat <- fm.lme$sigma
	S <- c(data.sp[1, 8, ])
	observed <- c(data.sp[3, 8, ]) == 1
	S.observed <- S[observed]
	ab.random <- fm.lme$coeff$random$subject
	a.random <- ab.random[order(as.numeric(row.names(ab.random))), 1:length(alpha.lme.int)]
	a.random <- t(t(a.random) + alpha.lme.int)
	b.random <- ab.random[order(as.numeric(row.names(ab.random))), (length(alpha.lme.int) + 1):(length(alpha.lme.int) + length(beta.lme.int))]
	b.random <- t(t(b.random) + beta.lme.int)
	a.random.observed <- a.random[observed, ]
	U.random.observed <- a.random.observed[, 1]
	lm.S.a <- lm(S.observed ~ U.random.observed) 
	gamma.int <- c(lm.S.a$coeff)
	sigma.s.hat <- sqrt(sum(lm.S.a$resid^2)/lm.S.a$df.residual)
    ## revise here when change model
	lm.b.a.0 <- lm(b.random[, 1] ~ 1)
	lm.b.a.1 <- lm(b.random[, 2] ~ -1 + a.random[, 2])
	beta.int.0 <- lm.b.a.0$coeff
	beta.int.1 <- lm.b.a.1$coeff
	beta.int <- c(beta.int.0, beta.int.1)
	if (trace) {
		cat("cov.a", cov.a.hat, "cov.b", cov.b.hat, "sigma.s", sigma.s.hat, "sigma.y", sigma.y.hat, "\n")
	}
    ## fit the shared-parameter model
	fm.sp <- NR.share.parameter(data = data, n.iter = n.iter, n.sample = n.sample, alpha.int = alpha.int, beta.int = beta.int, gamma.int = gamma.int, cov.a = cov.a.hat, cov.b = cov.b.hat, sigma.s = sigma.s.hat, sigma.y = sigma.y.hat, max.loglik = max.loglik, trace = trace)
	return(fm.sp)
}

#################################################
## fit shared-parameter model by fixed initial ##
#################################################
# set.seed(101)
fit.lme.sp.fixint <- function(data, n.iter, n.sample, max.loglik = T, trace = T) {
    ## revise here when change model
	alpha.int <- c(15.951, -1.983)
	beta.int <- c(-7.172, -0.23)
	gamma.int <- c(13.503, -0.403)
	cov.a.hat <- matrix(c(45.998, -3.727, -3.727, 1.379), 2, 2)
	cov.b.hat <- matrix(c(53.463, -8.115, -8.115, 5.872), 2, 2)
	sigma.y.hat <- 4.329
	sigma.s.hat <- 1.256
	if (trace) {
		cat("cov.a", cov.a.hat, "cov.b", cov.b.hat, "sigma.s", sigma.s.hat, "sigma.y", sigma.y.hat, "\n")
	}
    ## fit the shared-parameter model
	fm.sp <- NR.share.parameter(data = data, n.iter = n.iter, n.sample = n.sample, alpha.int = alpha.int, beta.int = beta.int, gamma.int = gamma.int,
           cov.a = cov.a.hat, cov.b = cov.b.hat, sigma.s = sigma.s.hat, sigma.y = sigma.y.hat, max.loglik = max.loglik, trace = trace)
	return(fm.sp)
}

###########################
## ENRICHD data analysis ##
###########################

enrichd.analysis <- function(data, min.mes, n.iter, n.sample, max.loglik = T, trace = T) {
    ## revise here when change model
	result <- array(NA, dim = c(n.iter + 1, 7, length(min.mes)))
	name.iter <- as.character(0:n.iter)
	alpha.item <- paste(rep("alpha", 2), as.character(0:1))
    ## revise here when change model
	beta.item <- paste(rep("beta", 2), c("00", "12"))
	gamma.item <- paste(rep("gamma", 2), as.character(0:1))
	name.item <- c("log-likelihood", alpha.item, beta.item, gamma.item)
	name.rc.sub <- paste(as.character(min.mes), rep("visits", length(min.mes)))
	dimnames(result) <- list(name.iter, name.item, min.mes)
	for (i in 1:length(min.mes)) {
		cat("at least", min.mes[i], "visits", "\n")
		data.sp <- data[, , data[2, 8, ] >= min.mes[i]]
		fit.sp <- try(fit.lme.sp(data = data.sp, n.iter = n.iter, n.sample = n.sample, max.loglik = max.loglik, trace = trace), silent = T)
		if (inherits(fit.sp, "try-error") == F) {
			result[, , i] <- fit.sp$whole
		}
	}
	return(result)
}

############################################
## ENRICHD data analysis by fixed initial ##
############################################

enrichd.analysis.fixint <- function(data, min.mes, n.iter, n.sample, max.loglik = T, trace = T) {
    ## revise here when change model
	result <- array(NA, dim = c(n.iter + 1, 7, length(min.mes)))
	name.iter <- as.character(0:n.iter)
	alpha.item <- paste(rep("alpha", 2), as.character(0:1))
    ## revise here when change model
	beta.item <- paste(rep("beta", 2), c("00", "12"))
	gamma.item <- paste(rep("gamma", 2), as.character(0:1))
	name.item <- c("log-likelihood", alpha.item, beta.item, gamma.item)
	name.rc.sub <- paste(as.character(min.mes), rep("visits", length(min.mes)))
	dimnames(result) <- list(name.iter, name.item, min.mes)
	for (i in 1:length(min.mes)) {
		cat("at least", min.mes[i], "visits", "\n")
		data.sp <- data[, , data[2, 8, ] >= min.mes[i]]
		fit.sp <- try(fit.lme.sp.fixint(data = data.sp, n.iter = n.iter, n.sample = n.sample, max.loglik = max.loglik, trace = trace), silent = T)
		if (inherits(fit.sp, "try-error") == F) {
			result[, , i] <- fit.sp$whole
		}
	}
	return(result)
}

##################
## bootstrap CI ##
##################

enrichd.bootstrap.CI <- function(data, min.mes, n.bts, n.iter, n.sample, max.loglik = T, trace = T) {
    ## revise here when change model
	result.bts <- array(NA, dim = c(1, 7, n.bts))
	data <- data[, , data[2, 8, ] >= min.mes]
	n.sub <- dim(data)[3]
	n.bts.eff <- 1
	cat("Bootstrap ")
	while (n.bts.eff <= n.bts) {
		sample.sub <- sample(1:n.sub, n.sub, replace = T)
		data.sp <- data[, , sample.sub]
		result.sp <- try(fit.lme.sp(data = data.sp, n.iter = n.iter, n.sample = n.sample, max.loglik = max.loglik, trace = trace), silent = T)
		if (inherits(result.sp, "try-error") == F) {
			cat(n.bts.eff)
			if (n.bts.eff %% 20 == 0) {
				cat("\n")
				if (n.bts.eff < n.bts) cat("Bootstrap ")
			}
			result.bts[, , n.bts.eff] <- result.sp$ultimate
			n.bts.eff <- n.bts.eff + 1
		}
	}
	cat("\n")
	return(result.bts)
}

###################################
## bootstrap CI by fixed initial ##
###################################

SPM.sub <- function(data=BDIdata, min.mes=5, n.boot=1000, n.iter=10, n.sample= 10000, max.loglik = T, trace = T) {
   
  medall.id<- 1:557
  medall.sp <- array(NA, dim = c(max(measurement), 8, length(medall.id)))
  
  for (i in 1:length(medall.id)) {
    ind.temp <- data[, 1] == medall.id[i]
    data <- data[ind.temp, ]                   # extract data from ith ID
    data <- matrix(data, ncol = dim(data)[2])
    n.mes <- dim(data)[1]                        # number of measurements
    t.point <- c(data[, 3])/30.5 # month         # vist time tij
    y <- c(data[, 2])                            # BDI at time t
    S <- data[1, 4]/30.5                         # med time Si
    R <- t.point - max(S, 0)                     # time after S if S<0, set 0
    lambda <- 1 * (t.point >= S)                 # med use indicator
    medall.sp[1, 8, i] <- S                      
    medall.sp[2, 8, i] <- n.mes
    medall.sp[3, 8, i] <- (S > min(t.point)) & (S <= max(t.point)) # with observed Si in [T1,Tn]
    medall.sp[4, 8, i] <- S == min(t.point)                        # used at first time T1. mostly the people start from baseline
    medall.sp[5, 8, i] <- S < min(t.point)                         # left censor
    medall.sp[6, 8, i] <- S > max(t.point)                         # right censor
    medall.sp[1:n.mes, 1, i] <- t.point                            # vector tij  of n.mes long
    medall.sp[1:n.mes, 2, i] <- y                                  # vector yij  of n.mes long
    medall.sp[1:n.mes, 3, i] <- 1                                  # n.mes of 1
    medall.sp[1:n.mes, 4, i] <- t.point                                  
    medall.sp[1:n.mes, 5, i] <- 1
    medall.sp[1:n.mes, 6, i] <- R                                  # vector Rij  of n.mes long
    medall.sp[1:n.mes, 7, i] <- lambda                             # vector delta_ij  of n.mes long
    
  }
  
 	result.bts <- array(NA, dim = c(1, 7, n.boot))
	data <- medall.sp[, , medall.sp[2, 8, ] >= min.mes]
	n.sub <- dim(data)[3]
	n.bts.eff <- 1
	cat("Bootstrap ")
	while (n.bts.eff <= n.boot) {
		sample.sub <- sample(1:n.sub, n.sub, replace = T)
		data.sp <- data[, , sample.sub]
		result.sp <- try(fit.lme.sp.fixint(data = data.sp, n.iter = n.iter, n.sample = n.sample, max.loglik = max.loglik, trace = trace), silent = T)
		if (inherits(result.sp, "try-error") == F) {
			cat(n.bts.eff)
			if (n.bts.eff %% 20 == 0) {
				cat("\n")
				if (n.bts.eff < n.boot) cat("Bootstrap ")
			}
			result.bts[, , n.bts.eff] <- result.sp$ultimate
			n.bts.eff <- n.bts.eff + 1
		}
	}
	cat("\n")
	return(result.bts)
}
