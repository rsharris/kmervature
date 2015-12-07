#===
#   kmervature-- model fitting to kmer abundance histograms
#   Copyright (C) 2015 Bob Harris
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#===
#
# The underlying functionality here was originally copied from KmerGenie,
# which can be found at http://kmergenie.bx.psu.edu.  This package extends
# the curve fitting used in KmerGenie, adding models for haploid+haploid and
# haploid+diploid mixtures and support for plotting results.
#
# Where functionality was borrowed from KmerGenie is indicated in comments.
# However, changes may have been made that deviate from specific KmerGenie
# functionality.
#
#===

#=== (not in KmerGenie) ===

read.coverage = function(hist.filename)
	{
	cov = read.table(hist.filename)[,1:2]
	cov[,2] = as.numeric(cov[,2])
	cov
	}


trim.coverage.tail = function(cov,copy.number.limit)
	{
	cov = cov[cov[,1] < copy.number.limit,]
	cov = cov[1:(max(which(cov[,2] != 0))), ]
	cov
	}


#=== (originally from KmerGenie's est-mean.r) ===

est.mean = function(d)
	{
	hc <- smooth(d[,2])

	# find first valley (or right before it)
	valley = hc[1]
	i = 2
	while (hc[i] < valley)
		{
		valley = hc[i]
		i = i + 1
		}

	# return max over the rest
	max.hist = hc[i]
	max.cov = i
	valley.pos = i
	while (i <= length(hc))
		{
		if (hc[i] > max.hist)
			{
			max.hist = hc[i]
			max.cov = i
			}
		i = i + 1
		}

	c(max.cov, valley.pos)
	}


#=== (originally from KmerGenie's est-params.r) ===

est.params <- function(d)
	{
	# First, find the valley and first coverage maximum using the smoothed histogram
	v <- est.mean(d)
	max.cov <- v[1]
	valley <- v[2]

	# Now find the median coverage value for everything past the valley
	vals <- d[,2]
	cvals <- cumsum(vals[(valley+1):length(vals)])
	before_median <- c(0,which(cvals < cvals[length(cvals)] / 2))
	mcov <- valley + max(before_median)

	# The new estimate of the mean coverage is the max of two estimates
	max.cov <- max(mcov, max.cov)

	# Try to estimate the deviation calculating the median absolute difference wrt
	# the calculated max.cov
	mvals <- c(vals[max.cov], vals[(max.cov+1):(2*max.cov-valley-1)] + rev(vals[(valley+1):(max.cov - 1)]))
	cvals <- cumsum(mvals)
	before_median <- c(0,which(cvals < cvals[length(cvals)] / 2))
	cov.sd <- 1.4826 * max(before_median - 1)
	cov.sd <- max(cov.sd,1.4826) # appears needed, sd=0 provokes crash in EMOpt

	# The initial estimate of erroroneous probability is just the ratio of k-mers
	# before the valley
	p.err <- sum(as.numeric(vals[1:valley])) / sum(as.numeric(vals))

	# FIXME: Estimate the skewness for skew normal
	# gamma <- min(0.9, abs(<sample skew>))
	# delta <- sqrt(pi / 2 * (gamma ^ (2/3)) / (gamma ^ (2/3) + ((4 - pi)/2) ^ (2/3)))
	# shape <- sign(<sample skew>) * delta / sqrt(1 - delta ^2)
	shape <- 0
	c(max.cov, cov.sd, p.err, shape)
	}


#=== (originally from KmerGenie's zeta.r) ===

zeta <- function (x)
	{
	a = 12
	k = 8
	B = c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
	ans = 0
	for (ii in 1:(a - 1)) ans = ans + 1/ii^x
	ans = ans + 1/((x - 1) * a^(x - 1)) + 1/(2 * a^x)
	term = (x/2)/a^(x + 1)
	ans = ans + term * B[1]
	for (mm in 2:k)
		{
		term = term * (x + 2 * mm - 2) * (x + 2 * mm - 3)/(a *
		a * 2 * mm * (2 * mm - 1))
		ans = ans + term * B[mm]
		}
	ans
	}

dzeta <- function(x, p, log. = FALSE)
	{
	x^(-p-1) / zeta(p + 1)
	}


#=== (originally from KmerGenie's model.r, which is for haploid) ===

perr <- function(x, shape, scale)
	{
	(1 + shape*(x-1)/scale) ^ (-1/shape) - (1 + shape*(x)/scale) ^ (-1/shape)
	}

truncate.copy <- 0

pgood <- function(x, max.copy, zp.copy, u.v, sd.v, shape.v)
	{
	# dcopy = heavy tail of the copy counts = heights of the zeta distribution with parameter zp.copy
	dcopy <- dzeta(1:max.copy, p=zp.copy)
	dcopy <- dcopy / sum(dcopy)

	# probs = matrix with max_abundance rows and max_copy_count columns, initialized to 0
	# it describes the mixture of distributions (after we sum all columns)

	if (truncate.copy == 0)
		{
		until.copy <- max.copy
		}
	else
		{
		until.copy <- truncate.copy
		}

	probs <- matrix(0, nrow = length(x), ncol=until.copy)

	for (copy in 1:until.copy)
		{
		# columns (copy) of = normal distribution correspond to copycount
		# dnorm(x,mean,sd) = height of the proba distribution at point x
		probs[,copy] = dcopy[copy] * 2 * dnorm(x, mean=(copy*u.v), sd=sqrt(copy)*sd.v) * pnorm((x - copy*u.v)*shape.v / (sqrt(copy)*sd.v))
		# FIXME: dnorm may be > 1 (consider the height of a normal distribution with extremly small sd), so this probabilty may exceed 1 with small datasets
		}

	rowSums(probs)
	}

model.haploid = function(params,cov,max.cov,max.copy)
	{
	zp.copy = params[1] # parameter of the zeta distribution, approximating the heavy tail
	p.e     = params[2] # "mixture parameter to determine which of the two [erroneous kmers and true kmers] distributions a k-mer coverage will be sampled from"
	shape.e = params[3] # pareto? shape
	scale.e = params[4] # pareto? scale
	u.v     = params[5] # abundance at copy_count=1
	sd.v    = params[6] # standard deviation of abundance at copy_count=1
	shape.v = params[7] # shape of the skewed normal at copy_count=1

	if (zp.copy <= 1 | shape.e < 0 | sd.v <= 0 | u.v <= 0 | p.e < 0 | p.e > 1 - 1e-9)
		return (list(like=-Inf))

	probs.err =  perr(cov[,1], shape.e, scale.e)
	probs.good = pgood(cov[,1], max.copy, zp.copy, u.v, sd.v, shape.v)

	kmers.probs <- p.e * probs.err + (1 - p.e) * probs.good
	logprobs <- log(kmers.probs)

	# for each abundance, multiply each log(sum(row)) by the number of kmers at that abundance
	return (list(like=-sum(logprobs*cov[,2]),
	             probs=kmers.probs,
	             probs.err=probs.err,
	             probs.good=probs.good))
	}

modelEM.haploid <- function(params,z,cov,max.cov,max.copy)
	{
	# note that p.e is replaced by z
	zp.copy = params[1]
	shape.e = params[2]
	scale.e = params[3]
	u.v     = params[4]
	sd.v    = params[5]
	shape.v = params[6]

	if (zp.copy < 1+1e-6 | shape.e < 1e-9 | sd.v < 1e-9 | u.v < 1e-9)
		return (list(like=-Inf))

	# error
	logprobs <- z * log(perr(cov[, 1], shape.e, scale.e))

	# no error
	logprobs <- logprobs + (1 - z)* log(pgood(cov[, 1], max.copy, zp.copy, u.v, sd.v, shape.v))

	-sum(logprobs * cov[, 2])
	}

EStep <- function(params,max.cov,max.copy)
	{
	zp.copy = params[1]
	p = params[2]
	shape.e = params[3]
	scale.e = params[4]
	u.v = params[5]
	sd.v = params[6]
	shape.v = params[7]

	pe <- p * perr(1:max.cov, shape.e, scale.e)
	pe / (pe + (1 - p) * pgood(1:max.cov, max.copy, zp.copy, u.v, sd.v, shape.v))
	}

#=== (originally from KmerGenie's model-diploid.r) ===

derr.old <- function(x, max.cov, shape.e)
	{
	derr <- diff(1 - 1/(1:(min(max(x), max.cov) + 1))^shape.e)
	derr[x]
	}

model.diploid = function(params,cov,max.cov,max.copy)
	{
	zp.copy = params[1] # parameter of the zeta distribution, approximating the heavy tail
	p.e     = params[2] # "mixture parameter to determine which of the two [erroneous kmers and true kmers] distributions a k-mer coverage will be sampled from"
	shape.e = params[3] # pareto? shape
	scale.e = params[4] # pareto? scale
	u.v     = params[5] # abundance at copy_count=1
	sd.v    = params[6] # standard dev of abundance at copy_count=1
	p.d     = params[7] # proportion of the genome which is heterozygous
	var.w   = params[8] # variance of abundance at copy_count=0.5; for some reason cannot fit with sd
	zp.copy.het = params[9] # heavy tail for het
	# TODO: use skewed normal here too

	if(zp.copy <= 1 | shape.e <= 0 | scale.e <= 0 | sd.v <= 0 | var.w <= 0 | zp.copy.het <= 1 | p.e < 1e-9 | p.e > 1 - 1e-9 | u.v <= 0)
		return(list(like=-Inf))

	#probs.err =  perr(cov[, 1], shape.e, scale.e) # for some reason that pareto function doesnt work with model-diploid
	probs.err =  derr.old(cov[, 1], max.cov, shape.e)

	shape.v = 0 # non-skewed
	probs.hom.good = pgood(cov[, 1], max.copy, zp.copy, u.v, sd.v, shape.v)
	probs.het.good = pgood(cov[, 1], max.copy, zp.copy.het, 0.5*u.v, 0.5*sqrt(var.w), shape.v)
	probs.good = p.d * probs.het.good + (1-p.d) * probs.hom.good
	kmers.probs <- p.e * probs.err + (1 - p.e) * probs.good
	logprobs <- log(kmers.probs)

	return (list(like=-sum(logprobs*cov[,2]),
	             probs=kmers.probs,
	             probs.err=probs.err,
	             probs.good=probs.good,
	             probs.het.good=probs.het.good,
	             probs.hom.good=probs.hom.good))
	}

#=== (not in KmerGenie) ===

model.enriched.hap.hap = function(params,cov,max.cov,max.copy)
	{
    model.enriched.hap.hap.params <<- params  # (global var) save 'em in case convergence fails

	zp.copy.y    = params[1]  # zeta shape for chrY part
	zp.copy.auto = params[2]  # zeta shape for autosomal part
	p.e          = params[3]  # mixture parameter, erroneous kmers vs true kmers
	shape.e      = params[4]  # pareto? shape (for erroneous kmers)
	scale.e      = params[5]  # pareto? scale (for erroneous kmers)
	p.y          = params[6]  # mixture parameter, chrY kmers vs autosomal kmers (relates to enrichment)
	u.y          = params[7]  # abundance at copy_count=1 for chrY part
	sd.y         = params[8]  # standard dev of abundance at copy_count=1 for chrY part
	shape.y      = params[9]  # shape of the skewed normal at copy_count=1 for chrY part
	u.auto       = params[10] # abundance at copy_count=1 for autosomal part
	sd.auto      = params[11] # standard dev of abundance at copy_count=1 for autosomal part

	eps = 1e-9
	if (zp.copy.y <= 1 | zp.copy.auto <= 1
	  | p.e < 0 | p.e > 1-eps | shape.e < 0 | scale.e <= 0
	  | p.y < eps | p.y > 1-eps
	  | u.y <= 0 | sd.y <= 0
	  | u.auto <= 0 | sd.auto <= 0)
		return (list(like=-Inf))

	p.good = 1 - p.e
	p.auto = 1 - p.y
	shape.auto = 0 # non-skewed

	probs.err       = perr(cov[,1], shape.e, scale.e)
	#probs.err      = derr.old(cov[,1], max.cov, shape.e)

	probs.y.good    = pgood(cov[,1], max.copy, zp.copy.y,    u.y,    sd.y,    shape.y)
	probs.auto.good = pgood(cov[,1], max.copy, zp.copy.auto, u.auto, sd.auto, shape.auto)

	probs.good  = p.y * probs.y.good + p.auto * probs.auto.good
	kmers.probs = p.e * probs.err    + p.good * probs.good
	logprobs <- log(kmers.probs)
    model.enriched.hap.hap.kmers.probs <<- kmers.probs  # (global var) save 'em in case convergence fails

	return (list(like=-sum(logprobs*cov[,2]),
	             probs=kmers.probs,
	             probs.err=probs.err,
	             probs.good=probs.good,
	             probs.y.good=probs.y.good,
	             probs.auto.good=probs.auto.good))
	}


model.enriched.hap.dip = function(params,cov,max.cov,max.copy)
	{
    model.enriched.hap.dip.params <<- params  # (global var) save 'em in case convergence fails

	zp.copy.y   = params[1]  # zeta shape for chrY part
	zp.copy.hom = params[2]  # zeta shape for homozygous autosome part
	zp.copy.het = params[3]  # zeta shape for heterozygous autosome part
	p.e         = params[4]  # mixture parameter, erroneous kmers vs true kmers
	shape.e     = params[5]  # pareto? shape (for erroneous kmers)
	scale.e     = params[6]  # pareto? scale (for erroneous kmers)
	p.y         = params[7]  # mixture parameter, chrY kmers vs autosomal kmers (relates to enrichment)
	u.y         = params[8]  # abundance at copy_count=1 for chrY part
	sd.y        = params[9]  # standard dev of abundance at copy_count=1 for chrY part
	shape.y     = params[10] # shape of the skewed normal at copy_count=1 for chrY part
	p.hom       = params[11] # mixture parameter, homozygous kmers vs heterozygous kmers (for autosomal part)
	u.hom       = params[12] # abundance at copy_count=1 for homozygous autosome part
	sd.hom      = params[13] # standard dev of abundance at copy_count=1 for homozygous autosome part
	var.het     = params[14] # variance of abundance at copy_count=0.5 for heterozygous autosome part

	eps = 1e-9
	if (zp.copy.y <= 1 | zp.copy.hom <= 1 | zp.copy.het <= 1
	  | p.e < 0 | p.e > 1-eps | shape.e < 0 | scale.e <= 0
	  | p.y < eps | p.y > 1-eps
	  | u.y <= 0 | sd.y <= 0
	  | p.hom < eps | p.hom > 1-eps
	  | u.hom <= 0 | sd.hom <= 0 | var.het <= 0)
		return (list(like=-Inf))

	p.good = 1 - p.e
	p.auto = 1 - p.y
	p.het  = 1 - p.hom
	shape.hom = 0 # non-skewed

	probs.err       = perr(cov[,1], shape.e, scale.e)
	#probs.err      = derr.old(cov[,1], max.cov, shape.e)

	probs.y.good    = pgood(cov[,1], max.copy, zp.copy.y,   u.y,       sd.y,              shape.y)
	probs.hom.good  = pgood(cov[,1], max.copy, zp.copy.hom, u.hom,     sd.hom,            shape.hom)
	probs.het.good  = pgood(cov[,1], max.copy, zp.copy.het, 0.5*u.hom, 0.5*sqrt(var.het), shape.hom)

	probs.auto.good = p.hom * probs.hom.good + p.het * probs.het.good
	probs.good      = p.y   * probs.y.good   + p.auto * probs.auto.good
	kmers.probs     = p.e   * probs.err      + p.good * probs.good
	logprobs <- log(kmers.probs)
    model.enriched.hap.dip.kmers.probs <<- kmers.probs  # (global var) save 'em in case convergence fails

	return (list(like=-sum(logprobs*cov[,2]),
	             probs=kmers.probs,
	             probs.err=probs.err,
	             probs.good=probs.good,
	             probs.y.good=probs.y.good,
	             probs.auto.good=probs.auto.good,
	             probs.het.good=probs.het.good,
	             probs.hom.good=probs.hom.good))
	}

#=== (originally from KmerGenie's fit-histogram.r) ===

# EM optimization

EMOpt <- function(cov,max.cov,max.copy,init, maxit = 1000, tol = 1e-9)
	{
	cpar <- init
	lik <- +Inf
	for (i in 1:maxit)
		{
		z <- EStep(cpar,max.cov,max.copy)

		p <- sum(z * cov[,2]) / sum(cov[,2])
		cpar = cpar[-2]
		#cat(cpar,"\n")
		opt <- optim(cpar,
		             modelEM.haploid,
		             z=z, cov=cov, max.cov=max.cov, max.copy=max.copy,
		             control=list(trace=0, maxit=10000))
		cpar <- append(opt$par, p, 1)
		clik <- model.haploid(cpar,cov,max.cov,max.copy)$like
		#cat("Log-likelihood:", clik, "\n")
		if (abs(clik - lik) < tol)
		break
		lik <- clik
		}

	cpar
	}

#=== (originally from KmerGenie's cutoff.r) ===

cutoffs <- function(probs.err, probs.good,n)
	{
	ratios = numeric(n)
	for (i in 1:n)
		{
		#ratios[i] = probs.err[i] / sum(probs.good[i])
		# suggested by claire: p.e*probs.err[i] / ((1-p.e)*all.probs.good[i])
		ratios[i] = (p.e*probs.err[i]) / ((1-p.e)*sum(probs.good[i]))
		cat(i," ",ratios[i],"\n")
		}

	return (ratios)
	}

#=== (not in KmerGenie) ===

kmer.histogram.fit <- function(model.type,cov,max.copy=3,params.init=NULL,
                               performFit=T,
                               zp.copy.est=NA,shape.e.est=NA,scale.e.est=NA,
                               p.y.est=NA,p.het.est=NA)
	{
	params = est.params(cov)
	cov.est   = params[1]
	sd.est    = params[2]
	p.e.est   = params[3]
	shape.est = params[4]

	if (is.na(zp.copy.est)) zp.copy.est = 3
	if (is.na(shape.e.est)) shape.e.est = 3
	if (is.na(scale.e.est)) scale.e.est = 1
	if (is.na(p.y.est))     p.y.est     = .9
	if (is.na(p.het.est))   p.het.est   = .2

	# filter extremes from kmers
	cov = trim.coverage.tail(cov,1.25*max.copy*cov.est)
	max.cov = max(cov[,1]) # maximum abundance in the histogram

	if (model.type == "haploid")
		{
		params.names <- c("zp.copy","p.e","shape.e","scale.e","u.v","sd.v","shape.v")
		if (is.null(params.init))
			{
			params.init = c(zp.copy.est,                        # zp.copy
			                p.e.est, shape.e.est, scale.e.est,  # p.e shape.e scale.e
			                cov.est, sd.est, shape.est)         # u.v sd.v shape.v
			}
		if (performFit)
			{
			kmer.histogram.fit.params.init <<- params.init  # (global)
			opt <- list(par = EMOpt(cov,max.cov,max.copy,params.init))
			params <- opt$par
			}
		}
	else if (model.type == "diploid")
		{
		params.names <- c("zp.copy","p.e","shape.e","scale.e","u.v","sd.v","p.d","var.w","zp.copy.het")
		if (is.null(params.init))
			{
			params.init = c(zp.copy.est,                        # zp.copy
			                p.e.est, shape.e.est, scale.e.est,  # p.e shape.e scale.e
			                cov.est, sd.est,                    # u.v sd.v
			                p.het.est, (sd.est)**2, 3)          # p.d var.w zp.copy.het
			}
		if (performFit)
			{
			kmer.histogram.fit.params.init <<- params.init  # (global)
			opt = optim(params.init,
			            function(x) model.diploid(x,cov,max.cov,max.copy)$like,
			            method="BFGS", control=list(trace=1, maxit=100), hessian=TRUE)
			params <- opt$par
			}
		}
	else if (model.type == "enriched-hap-hap")
		{
		params.names <- c("zp.copy.y","zp.copy.auto","p.e","shape.e","scale.e","p.y","u.y","sd.y","shape.y","u.auto","sd.auto")
		if (is.null(params.init))
			{
			p.auto.est = 1-p.y.est
			params.init = c(zp.copy.est, zp.copy.est,             # zp.copy.y zp.copy.auto
					        p.e.est, shape.e.est, scale.e.est,    # p.e shape.e scale.e
					        p.y.est, cov.est, sd.est, shape.est,  # p.y u.y sd.y shape.y
					        p.auto.est*cov.est, sd.est)           # u.auto sd.auto
			}
		if (performFit)
			{
			kmer.histogram.fit.params.init <<- params.init  # (global)
			opt = optim(params.init,
			            function(x) model.enriched.hap.hap(x,cov,max.cov,max.copy)$like,
			            method="BFGS", control=list(trace=1, maxit=100), hessian=TRUE)
			params <- opt$par
			}
		}
	else if (model.type == "enriched-hap-dip")
		{
		params.names <- c("zp.copy.y","zp.copy.hom","zp.copy.het","p.e","shape.e","scale.e","p.y","u.y","sd.y","shape.y","p.hom","u.hom","sd.hom","var.het")
		if (is.null(params.init))
			{
			p.auto.est = 1-p.y.est
			p.hom.est = 1-p.het.est
			u.hom.est = p.auto.est*p.hom.est*cov.est
			sd.hom.est = p.auto.est*p.hom.est*sd.est
			var.het.est = sd.hom.est**2
			params.init = c(zp.copy.est, zp.copy.est, zp.copy.est,          # zp.copy.y zp.copy.hom zp.copy.het
					        p.e.est, shape.e.est, scale.e.est,              # p.e shape.e scale.e
					        p.y.est, cov.est, sd.est, shape.est,            # p.y u.y sd.y shape.y
					        p.hom.est, u.hom.est, sd.hom.est, var.het.est)  # p.hom u.hom sd.hom var.het
			}
		if (performFit)
			{
			kmer.histogram.fit.params.init <<- params.init  # (global)
			opt = optim(params.init,
			            function(x) model.enriched.hap.dip(x,cov,max.cov,max.copy)$like,
			            method="BFGS",
						control=list(trace=1, maxit=100), hessian=TRUE)
			params <- opt$par
			}
		}
	else
		{
		print (paste("unknown model.type:",model.type))
		return (NULL)
		}

	if (!performFit)
		{
		return (list(cov = cov,
		             max.cov = max.cov,
		             max.copy = max.copy,
		             params.init = params.init,
		             params.names = params.names,
		             opt = NA))
		}

	list(cov = cov,
	     max.cov = max.cov,
	     max.copy = max.copy,
	     params.init = params.init,
	     params.names = params.names,
	     params = opt$par,
	     opt = opt)
	}


model.haploid.plot <- function(fit,main)
	{
	cov = fit$cov
	mm = model.haploid(fit$params,cov,fit$max.cov,fit$max.copy)
	p.err  = fit$params[2]
	p.good = 1 - p.err
	numKmers = sum(cov[,2])

	legText = c("actual","error","good")
	legCol  = c("black","red","blue")
	legLty  = c(2,1,1)
	legLwd  = c(4,1,4)

	plot  (cov[,1],p.err*mm$probs.err*numKmers,col=legCol[2],lty=legLty[2],lwd=legLwd[2],
		   type="l",log="y",ylim=c(1,max(cov[,2])),
		   main=main,
		   xlab="number of times a distinct kmer occurs",
		   ylab="number of kmers (log)")
	lines (cov[,1],p.good*mm$probs.good*numKmers,col=legCol[3],lty=legLty[3],lwd=legLwd[3])
	lines (cov,col=legCol[1],lty=legLty[1],lwd=legLwd[1])

	legend("bottomright",legText,lty=legLty,lwd=legLwd,col=legCol)
	}


model.diploid.plot <- function(fit,main)
	{
	cov = fit$cov
	mm = model.diploid(fit$params,cov,fit$max.cov,fit$max.copy)
	p.err  = fit$params[2]
	p.good = 1 - p.err
	p.het  = fit$params[7]
	p.hom  = 1 - p.het
	numKmers = sum(cov[,2])

	legText = c("actual","error","good","homozygous","heterozygous")
	legCol  = c("black","red","blue","#FF8000","#008000")
	legLty  = c(2,1,1,2,2)
	legLwd  = c(4,1,4,1,1)

	plot  (cov[,1],p.err*mm$probs.err*numKmers,col=legCol[2],lty=legLty[2],lwd=legLwd[2],
		   type="l",log="y",ylim=c(1,max(cov[,2])),
		   main=main,
		   xlab="number of times a distinct kmer occurs",
		   ylab="number of kmers (log)")
	lines (cov[,1],p.good*mm$probs.good*numKmers,col=legCol[3],lty=legLty[3],lwd=legLwd[3])
	lines (cov[,1],p.good*p.hom*mm$probs.hom.good*numKmers,col=legCol[4],lty=legLty[4],lwd=legLwd[4])
	lines (cov[,1],p.good*p.het*mm$probs.het.good*numKmers,col=legCol[5],lty=legLty[5],lwd=legLwd[5])
	lines (cov,col=legCol[1],lty=legLty[1],lwd=legLwd[1])

	legend("bottomright",legText,lty=legLty,lwd=legLwd,col=legCol)
	}


model.enriched.hap.hap.plot <- function(fit,main)
	{
	cov = fit$cov
	mm = model.enriched.hap.hap(fit$params,cov,fit$max.cov,fit$max.copy)
	p.err  = fit$params[3]
	p.good = 1 - p.err
	p.y    = fit$params[6]
	p.auto = 1 - p.y
	numKmers = sum(cov[,2])

	legText = c("actual","error","non-error","haploid (chrY)","haploid (chr8)")
	legCol  = c("black","red","blue","green","red")
	legLty  = c(2,1,1,2,1)
	legLwd  = c(4,2,4,2,2)

	plot  (cov[,1],p.err*mm$probs.err*numKmers,col=legCol[2],lty=legLty[2],lwd=legLwd[2],               # actual
		   type="l",log="y",ylim=c(1,max(cov[,2])),
		   main=main,
		   xlab="number of times a distinct kmer occurs",
		   ylab="number of kmers (log)")
	lines (cov[,1],p.good       *mm$probs.good     *numKmers,col=legCol[3],lty=legLty[3],lwd=legLwd[3]) # non-error
	lines (cov[,1],p.good*p.y   *mm$probs.y.good   *numKmers,col=legCol[4],lty=legLty[4],lwd=legLwd[4]) # haploid (chrY)
	lines (cov[,1],p.good*p.auto*mm$probs.auto.good*numKmers,col=legCol[5],lty=legLty[5],lwd=legLwd[5]) # haploid (chr8)
	lines (cov,col=legCol[1],lty=legLty[1],lwd=legLwd[1])                                               # error

	legend("bottomright",legText,lty=legLty,lwd=legLwd,col=legCol)
	}


model.enriched.hap.dip.plot <- function(fit,main)
	{
	cov = fit$cov
	mm = model.enriched.hap.dip(fit$params,cov,fit$max.cov,fit$max.copy)
	p.err  = fit$params[4]
	p.good = 1 - p.err
	p.y    = fit$params[7]
	p.auto = 1 - p.y
	p.hom  = fit$params[11]
	p.het  = 1 - p.hom
	numKmers = sum(cov[,2])

	legText = c("actual","error","non-error","haploid (chrY)","diploid (chr8)","diploid homozygous","diploid heterozygous")
	legCol  = c("black","red","blue","green","red","#FF8000","#008000")
	legLty  = c(2,1,1,2,2,2,2)
	legLwd  = c(4,2,4,2,2,1,1)

	plot  (cov[,1],p.err*mm$probs.err*numKmers,col=legCol[2],lty=legLty[2],lwd=legLwd[2],                     # actual
		   type="l",log="y",ylim=c(1,max(cov[,2])),
		   main=main,
		   xlab="number of times a distinct kmer occurs",
		   ylab="number of kmers (log)")
	lines (cov[,1],p.good             *mm$probs.good     *numKmers,col=legCol[3],lty=legLty[3],lwd=legLwd[3]) # non-error
	lines (cov[,1],p.good*p.y         *mm$probs.y.good   *numKmers,col=legCol[4],lty=legLty[4],lwd=legLwd[4]) # haploid (chrY)
	lines (cov[,1],p.good*p.auto      *mm$probs.auto.good*numKmers,col=legCol[5],lty=legLty[5],lwd=legLwd[5]) # diploid (chr8)
	lines (cov[,1],p.good*p.auto*p.hom*mm$probs.hom.good *numKmers,col=legCol[6],lty=legLty[6],lwd=legLwd[6]) # diploid homozygous
	lines (cov[,1],p.good*p.auto*p.het*mm$probs.het.good *numKmers,col=legCol[7],lty=legLty[7],lwd=legLwd[7]) # diploid heterozygous
	lines (cov,col=legCol[1],lty=legLty[1],lwd=legLwd[1])                                                     # error

	legend("bottomright",legText,lty=legLty,lwd=legLwd,col=legCol)
	}


kmer.histogram.plot <- function(cov,main,xlab=NA,ylab=NA,
                                col="black",lty=1,lwd=1)
	{
	if (is.na(xlab))
		{
		xlab = "number of times a distinct kmer occurs"
		}

	if (is.na(ylab))
		{
		ylab = "number of kmers (log)"
		}

	plot(cov,type="l",log="y",
	     main=main,xlab=xlab,ylab=ylab,
	     col=col,lty=lty,lwd=lwd)
	}
