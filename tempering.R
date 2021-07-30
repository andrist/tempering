#!/usr/bin/Rscript

# required libraries and helper functions (common.R)
path = "."
if (i<-grep("--file=", commandArgs())) {
	path = dirname(sub("--file=", "", commandArgs()[i]))
}
source(paste(path, '/', 'common.R', sep=""), chdir=T)

# list of command line options for the user
parser <- OptionParser(usage="        %prog [options]

examples:
        %prog --range 0.5,2.8 -n 16           # geometric series from 0.5 to 2.8
        %prog data.txt --estimator energy     # use energy estimator on data.txt
        %prog data.txt --combined 5,4,0,1,0   # combine multiple weighted methods",
	option_list = list(
	make_option(c('-e', '--estimator'), action='store',
		help='estimator to use (flow,energy,accrate,geometric,linear)'),
	make_option(c('-c', '--combined'), action='store',
		help='weighted combination of estimators (order as above, default=5,4,0,1,0)'),
	make_option(c('-l', '--labels'), action='store', default='pts,ptu,ptd',
		help='header used for acceptance, up-flow and down-flow'),
	make_option(c('-r', '--range'), action='store', default=NULL,
		help='temperature range to cover: Tmin,Tmax'),
	make_option(c('-v', '--verbose'), action='store_true', default=F,
		help='show verbose output on optimization process'),
	make_option(c('-n', '--size'), action='store', default=NULL,
		help='target number of temperatures in the set')
))
config <- parse_args(parser, positional_arguments=TRUE)
opt <- config$options
arg <- config$args

if (length(arg) == 0 && (is.null(opt$range) || is.null(opt$size))) {
	print_help(parser)
	quit()
}

# script configuration
estim <- c('flow', 'energy', 'accrate', 'geometric', 'linear')
color <- c('blue', 'red', 'dark green', 'dark gray', 'gray')
weight <- c(5, 4, 0, 1, 0)
if (!is.null(opt$estimator)) {
	found=F
	for (i in 1:5) {
		if (estim[i] == opt$estimator) {
			weight=rep(0, 5)
			weight[i]=1
			found=T
		}
	}
	if (!found) {
		stop(paste("estimator '", opt$estimator, "' does not exist.", sep=""))
	}
}
if (!is.null(opt$combined)) {
	weight <- as.numeric(strsplit(opt$combined, ",")[[1]])
	weight <- c(weight, rep(5, 0))
	weight <- weight[1:5]
}

# handle dat input from file
if (length(arg)==0) {
	dat <- NULL
  weight[1] = weight[2] = weight[3] = 0;
} else if (length(arg)==1) {
	dat <- read.table(arg[1], header=T)
	if (!'ptu' %in% colnames(dat)) {
		if ('ptd' %in% colnames(dat)) {
			dat$ptu = 1 - dat$ptd
		} else {
			dat$ptu = (0:(-(nrow(dat)-1)))/(nrow(dat)-1)
			if (opt$verbose || !is.null(opt$weight) && weight[1] > 0) {
				write("column ptu or ptd not found, flow estimator", stderr())
			}
			weight[1] = 0
		}
	}
	if (!'ptd' %in% colnames(dat)) {
		dat$ptd = 1-dat$ptu
	}
	if (!'e' %in% colnames(dat)) {
		dat$e = (0:(-(nrow(dat)-1)))/(nrow(dat)-1)
		if (opt$verbose || !is.null(opt$weight) && weight[2] > 0) {
			write("column e not found, energy estimator turned off", stderr())
		}
		weight[2] = 0
	}
	if (!'pts' %in% colnames(dat)) {
		dat$pts=rep(1,nrow(dat))
		if (opt$verbose || !is.null(opt$weight) && weight[3] > 0) {
			write("column pts not found, acceptance estimator turned off", stderr())
		}
		weight[3] = 0
	}
} else {
	stop('too many arguments, specify at most 1 dat file')
}
if (sum(weight)<=0) {
	write("total weight is zero, defaulting to geometric series", stderr())
	weight = c(0, 0, 0, 1, 0)
}
print(weight)

# command line options
if (!is.null(dat)) {
	Trng <- range(dat$T)
}
if (!is.null(opt$range)) {
	Trng <- as.numeric(strsplit(opt$range, ",")[[1]])
}
if (is.null(Trng)) {
	stop('error: T range is null: no file or --range option specified')
}
Tmin <- min(Trng)
Tmax <- max(Trng)
N <- length(dat$T)
if (!is.null(opt$size)) {
	N <- as.numeric(opt$size)
}
Tset <- rep(Tmin,N)
Tset[N] <- Tmax

# densities for the various methods
n <- length(dat$T)

if (weight[1]>0) {
	# flow estimator
	sm_eta <- etas <- diff(dat$ptd/(dat$ptu+dat$ptd))/diff(dat$T)
	for (threshold in 0.9^(20:100)) {
		mid_eta <- mids(sm_eta);
		if (sm_eta[1]<threshold)
			sm_eta[1] <- mean(c(sm_eta[1], mid_eta[1]))
		for (i in 4:n-2) if (sm_eta[i]<threshold)
			sm_eta[i] = mean(mid_eta[(i-1):i])
		if (sm_eta[n-1]<threshold)
			sm_eta[n-1] <- mean(mid_eta[n-2], sm_eta[n-1])
	}
	eta <- approxfun(mids(dat$T), sm_eta, rule=2)
} else {
	etas <- diff(dat$T)*0;
	eta <- function(x) { return(x*0); }
}

# acceptance estimator (to be used iteratively)
if (weight[3]>0) {
	accr <- approxfun(mids(dat$T), (1-dat[['pts']][2:n-1])/diff(dat$T), rule=2)
} else {
	accr <- function(x) { return(x*0); }
}

# energy estimator (my own estimator)
if (weight[2] > 0) {
	eN <- max(50, N)
	es <- approxfun(dat$T, dat$e, rule=2);
	eTset <- function(step_) {
		cost <- function(T_) (es(T)-es(T_))*(1./T-1./T_)-step_
		Ts <- rep(0,eN); Ts[1] <- T <- Tmin
		for (i in 2:eN) {
			T2 <- T; while( cost(T2)>0 && T2<Tmax*2) T2 <- T2*2;
			if (T2>Tmax*2) return(c(Tmax*2-step_))
			Ts[i] <- T <- uniroot(cost, c(T, T2), tol=1e-10)$root
		}
		return(Ts)
	}
	eTmax <- function(step) Tmax-max(eTset(step))
	step <- uniroot(eTmax, c(-1e-1, -1e-15), tol=1e-15)$root
	eTs <- eTset(step)
	edens <- approxfun(mids(eTs), 1./diff(eTs)/(eN-1), rule=2)
} else {
	edens <- function(x) { return(x*0); }
}

# list of estimators
estim <- list(
	function(T) eta(T),
	function(T) edens(T),
	function(T) accr(T),
	function(T) 1/T,
	function(T) rep(1, length(T))
)

weight <- weight/sum(weight)
normal <- rep(1., length(estim))
for (i in 1:length(estim)) {
	normal[i] <- 1./integrate(estim[[i]], Tmin, Tmax, rel.tol=1e-3)$value
}

combined <- function(T) {
	dens <- rep(0,length(T))
	for (i in 1:length(estim)) {
    if (weight[i]>0) {
		  dens = dens + estim[[i]](T)*normal[i]*weight[i]
    }
  }
	return(dens)
}

# reduce high frequency noise by resampling
sn = 50
sT <- 0:sn/sn*(Tmax - Tmin) + Tmin
sD <- combined(sT)
combined <- splinefun(sT, sD)

# generate plot comparing the different densities
Tk = 0:500/500*(Tmax - Tmin) + Tmin
rng <- range(c(eta(Tk), edens(Tk), combined(Tk)))
plot(Tk, combined(Tk), col='black', type='l', ylim=rng,
	 main='Temperature Set Density')
for (i in 1:length(estim)) {
  if (weight[i] > 0) {
    lines(Tk, normal[i]*estim[[i]](Tk), col=color[i], lty=2)
  }
}
points(mids(dat$T), etas)
lines(mids(dat$T), 1/diff(dat$T)/(n-1), lty=3)
legend(Tmin+(Tmax-Tmin)*2/3, rng[2],
	   c('combined', paste(estim, sprintf("(%.2f)", weight)), 'previous'),
	   col=c('black', color, 'black'),
	   lty=c(1, rep(2, 5), 3))

intc <- integrate(combined, Tmin, Tmax, rel.tol=1e-3)$value
# generate and write temperature set to stdout
for (i in 3:N-1) {
  f <- function(Te) integrate(combined, Tset[i-1], Te, rel.tol=1e-3)$value-intc/(N - 1)
	Tset[i] <- uniroot(f, c(Tset[i-1], Tset[N]), tol=1e-6)$root
}
sTset <- paste(sprintf("%.4f", Tset-4.9e-5), collapse=" ")
write(paste("T=", sTset, sep=""), stdout())

# generate plot comparing old/new temperature set
if (F && !is.null(dat)) {
	rng <- range(diff(Tset)/mean(diff(Tset)));
	rng <- range(c(rng, diff(dat$T)/mean(diff(dat$T))))
	plot(mids(Tset), diff(Tset)/mean(diff(Tset)), ylim=rng,
		   main='Spacing for the new Temperature Set')
	lines(mids(Tset), diff(Tset)/mean(diff(Tset)))
	points(mids(dat$T), diff(dat$T)/mean(diff(dat$T)), col='red')
	lines(mids(dat$T), diff(dat$T)/mean(diff(dat$T)), col='red')
	legend(Tmin+(Tmax-Tmin)*2/3, .5, c('suggested', 'previous'),
		   col=c('black', 'red'), lty=c(1, 1))
}

