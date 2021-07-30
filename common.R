#!/usr/bin/Rscript

#Note: using UTF-8 causes problem with X11 on mac
invisible(Sys.setlocale("LC_CTYPE", "C"));

#These scripts require two packages, search & install
install.dependencies <- function(dependencies) {
	for (dep in dependencies) {
		r <- 'http://cran.R-project.org'
		if (!is.element(dep, installed.packages())) {
			write(paste('warning: required package', dep,
						      'not found, installing...'), stderr())
			invisible(install.packages(dep, repos=r))
			if (is.element(dep, installed.packages())) {
				write("installation successful", stderr())
			} else {
        stop(paste('installing', dep, 'failed,',
							     'please install manually'))
      }
		}
	}
}
libs <- c('getopt', 'optparse', 'boot')
install.dependencies(libs)
suppressPackageStartupMessages( library(optparse) )
# suppressPackageStartupMessages( library(boot) )
library(boot)

# little helper functions for rearranging data
mids <- function(x) x[2:length(x)] - diff(x)/2
cut_ends <- function(x, n=1) {
	n = max(0, n)
	if (length(x) <= 2*n) {
    return(x[c()])
  }
	return (x[(2*n + 1):length(x) - n])
}
add_ends <- function(x) c(x[1], x, x[length(x) - 1])
lowpass <- function(x, n=1) {
	for (i in 1:n) {
    x<-cut_ends(filter(add_ends(x), c(.25, .5, .25), side=2))
  }
	return (x)
}
interleave <- function(x, y){
	if (is.data.frame(x) || is.data.frame(y)) {
		return (cbind(x, y)[order(c(1:length(colnames(x)), 1:length(colnames(y))))])
	}else{
		return (c(x, y)[order(c(1:length(x), 1:length(y)))])
	}
}

# find directories
list.data.dirs <- function(path="./", pattern=".*\\.data$"){
	files <- list.files(path=path, pattern=pattern, recursive=T, full.names=T)
	return (unique(dirname(files)))
}
read.data.files <- function(path="./", pattern=".*\\.data$",
							              select.col='te', select='max', select.global=T) {
	files <- list.files(path=path, pattern=pattern, recursive=F, full.names=T)
	first <- T
	for (f in files) {
		d <- read.table(f, header=T)
		if (select.col %in% colnames(d)) {
			if (select=='max') {
				d <- d[d[[select.col]] == max(d[[select.col]]),]
				if (select.global && !first && 
				    max(d[[select.col]]) > max(dat[[select.col]])) first=T;
			} else if (select=='min') {
				d <- d[d[[select.col]] == min(d[[select.col]]),]
				if (select.global && !first && 
					  min(d[[select.col]]) < min(dat[[select.col]])) first=T;
			} else {
				d <- d[d[[select.col]] == select];
			}
		}
		if (first) {
      dat <- d
    } else {
      dat <- rbind(dat, d)
    }
		first <- F
	}
	return(dat)
}
