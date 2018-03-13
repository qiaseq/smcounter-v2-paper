# re-estimate background error from duplex-seq data 
# piece-wise density estimation with Pareto tail and beta distribution
# Chang Xu, 13NOV2017

rm(list=ls())
setwd('C:/Users/xuc/Documents/duplexTag/error_reestimated')
library(fitdistrplus)
library(Hmisc)
library(gPdtest)
library(dplyr)
library(laeken)
library(rmutil)

options(stringsAsFactors=F)
set.seed(732017)
nsim <- 50000
select <- dplyr::select
filter <- dplyr::filter

# function to find reverse nucleotide
rev <- function(x){
  if(x=='a') y <- 't'
  else if (x=='c') y <- 'g'
  else if (x=='g') y <- 'c'
  else if (x=='t') y <- 'a'
  else if (x=='A') y <- 'T'
  else if (x=='C') y <- 'G'
  else if (x=='G') y <- 'C'
  else if (x=='T') y <- 'A'
  else y <- 'n'  
  return(y)
}

min.mtDepth <- 1000
pMax <- 0.01
# this file contains HC region (v3.3.2), homozygous reference, filter passed, non-repetitive region sites
bkg1 <- read.delim('intermediate/bkg.duplex_oneReadPairUMIsDropped.txt', header=T)
bkg2 <- read.delim('intermediate/bkg.duplex_oneReadPairUMIsIncluded.txt', header=T)

types <- c('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')

#####################################################
# excluding 1 read MTs
#####################################################
dat <- bkg1
for(type in types){
  ref <- unlist(strsplit(type, split=''))[1]
  alt <- unlist(strsplit(type, split=''))[2]  
  tmp <- filter(dat, (REF==ref & negStrand > min.mtDepth) | (REF==rev(ref) & posStrand > min.mtDepth)) %>% 
    mutate(counts = eval(parse(text=paste0(ref, '.', alt))),
           all = ifelse(REF==ref, negStrand, posStrand),
           p = counts / all ) %>%
    filter(p < pMax)
  assign(paste0('p', ref, alt), tmp$p)
}

# 8 lower-error subs
tmp <- mutate(dat, sum8 = A.C + C.A + A.T + T.A + C.G + G.C + G.T + T.G, total8 = AllSMT * 2, p8 = sum8 / total8) %>% 
  filter(total8 > 80000 & p8 < .001)
plot(density(tmp$p8))
summary(tmp$p8)

# # Beta distribution + P(0)
# fit.beta.mle.ag <- fitdist(pAG[pAG>0], distr='beta', method='mle')
# fit.beta.mle.ga <- fitdist(pGA[pGA>0], distr='beta', method='mle')
# fit.beta.mle.ct <- fitdist(pCT[pCT>0], distr='beta', method='mle')
# fit.beta.mle.tc <- fitdist(pTC[pTC>0], distr='beta', method='mle')
# 
# p0 <- c(length(pAG[pAG==0])/length(pAG), length(pGA[pGA==0])/length(pGA), length(pCT[pCT==0])/length(pCT), length(pTC[pTC==0])/length(pTC))
# shape1 <- c(fit.beta.mle.ag$estimate[1], fit.beta.mle.ga$estimate[1], fit.beta.mle.ct$estimate[1], fit.beta.mle.tc$estimate[1])
# shape2 <- c(fit.beta.mle.ag$estimate[2], fit.beta.mle.ga$estimate[2], fit.beta.mle.ct$estimate[2], fit.beta.mle.tc$estimate[2])
# type <- c('A/G', 'G/A', 'C/T', 'T/C')
# top4 <- data.frame(type, p0, shape1, shape2)
# bkg.error <- list(top4.exclude.1rpUMI=top4)
# 
# # generate random numbers according to the fitted distributions
# r.beta.ag <- c(rep(0, round(nsim*p0[1])), rbeta(nsim-round(nsim*p0[1]), shape1[1], shape2[1]))
# r.beta.ga <- c(rep(0, round(nsim*p0[2])), rbeta(nsim-round(nsim*p0[2]), shape1[2], shape2[2]))
# r.beta.ct <- c(rep(0, round(nsim*p0[3])), rbeta(nsim-round(nsim*p0[3]), shape1[3], shape2[3]))
# r.beta.tc <- c(rep(0, round(nsim*p0[4])), rbeta(nsim-round(nsim*p0[4]), shape1[4], shape2[4]))
# 
# # evaluate the fit -- QQplot
# png('QQplot.beta.exclude_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   r <- eval(parse(text=paste0('r.beta.', tolower(type)))) 
#   p <- eval(parse(text=paste0('p', type)))
#   qqplot(p, r, xlab='observed', ylab='fitted', main=type)
#   abline(a=0, b=1, col='red', lty=2)
# }
# dev.off()
# 
# # evaluate the fit -- density function in the tail
# png('tail.beta.exclude_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   p <- eval(parse(text=paste0('p', type)))
#   plot(density(p[p>0]), xlim=c(.0008, .002), ylim=c(0, 100), main=type, xlab='')
#   fit <- eval(parse(text=paste0('fit.beta.mle.', tolower(type))))
#   a <- fit$estimate[1]
#   b <- fit$estimate[2]
#   curve(dbeta(x, a, b), add=T, col='red', lty=2)
# }
# dev.off()
# 
# # evaluate the fit -- ecdf
# png('ECDF.beta.exclude_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   r <- eval(parse(text=paste0('r.beta.', tolower(type)))) 
#   p <- eval(parse(text=paste0('p', type)))
#   plot(ecdf(p), xlab='', main=type, ylim=c(0.95, 1))
#   lines(ecdf(r), lty=2, col='red')
#   legend(x='bottomright', legend=c('observed', 'fitted'), lty=1:2, col=c('black', 'red'))
# }
# dev.off()


#######################  Beta distribution + Pareto tail  ######################
types4 <- c('AG', 'GA', 'CT', 'TC')

x0.AG <- 1.2e-3
x0.GA <- 1.2e-3
x0.CT <- 5e-4
x0.TC <- 5e-4

ul.AG <- 0.0015
ul.GA <- 0.0025
ul.CT <- 0.0020
ul.TC <- 1

# fit beta distribution; 0's are imputed
x0s <- thetas <- shape1s <- shape2s <- pBetas <- rep(NA, 4)

for(i in 1:4){
  type <- types4[i]
  # fit pareto tails
  ul <- eval(parse(text=paste0('ul.', type)))
  p.tmp <- eval(parse(text=paste0('p', type)))
  x0 <- eval(parse(text=paste0('x0.', type)))
  pareto.theta <- paretoTail(p.tmp[p.tmp < ul], x0 = x0)$theta
  # fit beta distribution
  p.min <- min(p.tmp[p.tmp > 0])
  p.tmp1 <- ifelse(p.tmp == 0, runif(1, 0, p.min), p.tmp)
  fit <- fitdist(p.tmp1, distr='beta', method='mle')
  shape1 <- fit$estimate[1]
  shape2 <- fit$estimate[2]
  
  # simulate beta's
  tmp <- rbeta(5*nsim, shape1, shape2)
  tmp <- tmp[tmp <= x0]
  p.beta <- sum(p.tmp <= x0)/length(p.tmp)
  r.beta <- sample(tmp, round(nsim*p.beta))
  
  # simulate pareto's; remove extremely large data points
  r.pare <- x0 * runif(round(nsim * (1-p.beta)))^(-1/pareto.theta)
  p.max <- max(p.tmp)
  r.pare <- ifelse(r.pare <= ul, r.pare, runif(1, 0, p.max))
  assign(paste0('r.', tolower(type)), c(r.beta, r.pare)) 
  
  # save parameters
  x0s[i] <- x0
  thetas[i] <- pareto.theta
  shape1s[i] <- shape1
  shape2s[i] <- shape2
  pBetas[i] <- p.beta
}

parameters.exclude_1rpUMI <- data.frame(types4, x0s, thetas, shape1s, shape2s, pBetas)
colnames(parameters.exclude_1rpUMI) <- c('type', 'x0', 'theta', 'shape1', 'shape2', 'pBeta')
r.exc1.ag <- r.ag
r.exc1.ga <- r.ga
r.exc1.ct <- r.ct
r.exc1.tc <- r.tc

# evaluate the fit -- QQplot
# png('QQplot.piecewise.exclude_1rpUMI.png', height=1200, width=1200, res=200)
png('tmp.QQplot.piecewise.exclude_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  r <- eval(parse(text=paste0('r.', tolower(type)))) 
  p <- eval(parse(text=paste0('p', type)))
  qqplot(p, r, xlab='observed', ylab='fitted', main=type)
  abline(a=0, b=1, col='red', lty=2)
}
dev.off()

# evaluate the fit -- full density
png('density.piecewise.exclude_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  p <- eval(parse(text=paste0('p', type)))
  plot(density(p), xlim=c(0, .002), ylim=c(0, 5000), main=type, xlab='')
  
  fit <- eval(parse(text=paste0('r.', tolower(type))))
  lines(density(fit), col='red', lty=2)
}
dev.off()

# evaluate the fit -- density function in the tail
png('tail.piecewise.exclude_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  p <- eval(parse(text=paste0('p', type)))
  plot(density(p), xlim=c(0.0008, .002), ylim=c(0, 100), main=type, xlab='')
  
  fit <- eval(parse(text=paste0('r.', tolower(type))))
  lines(density(fit), col='red', lty=2)
}
dev.off()

# evaluate the fit -- ecdf
png('ECDF.piecewise.exclude_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  r <- eval(parse(text=paste0('r.', tolower(type)))) 
  p <- eval(parse(text=paste0('p', type)))
  plot(ecdf(p), xlab='', main=type, ylim=c(0.95, 1))
  lines(ecdf(r), lty=2, col='red')
  legend(x='bottomright', legend=c('observed', 'fitted'), lty=1:2, col=c('black', 'red'))
}
dev.off()



#####################################################
# including 1 read MTs
#####################################################
dat <- bkg2
for(type in types){
  ref <- unlist(strsplit(type, split=''))[1]
  alt <- unlist(strsplit(type, split=''))[2]  
  tmp <- filter(dat, (REF==ref & negStrand > min.mtDepth) | (REF==rev(ref) & posStrand > min.mtDepth)) %>% 
    mutate(counts = eval(parse(text=paste0(ref, '.', alt))),
           all = ifelse(REF==ref, negStrand, posStrand),
           p = counts / all ) %>%
    filter(p < pMax)
  assign(paste0('p', ref, alt), tmp$p)
}

# # Beta distribution + P(0)
# fit.beta.mle.ag <- fitdist(pAG[pAG>0], distr='beta', method='mle')
# fit.beta.mle.ga <- fitdist(pGA[pGA>0], distr='beta', method='mle')
# fit.beta.mle.ct <- fitdist(pCT[pCT>0], distr='beta', method='mle')
# fit.beta.mle.tc <- fitdist(pTC[pTC>0], distr='beta', method='mle')
# 
# p0 <- c(length(pAG[pAG==0])/length(pAG), length(pGA[pGA==0])/length(pGA), length(pCT[pCT==0])/length(pCT), length(pTC[pTC==0])/length(pTC))
# shape1 <- c(fit.beta.mle.ag$estimate[1], fit.beta.mle.ga$estimate[1], fit.beta.mle.ct$estimate[1], fit.beta.mle.tc$estimate[1])
# shape2 <- c(fit.beta.mle.ag$estimate[2], fit.beta.mle.ga$estimate[2], fit.beta.mle.ct$estimate[2], fit.beta.mle.tc$estimate[2])
# type <- c('A/G', 'G/A', 'C/T', 'T/C')
# top4 <- data.frame(type, p0, shape1, shape2)
# bkg.error <- list(top4.exclude.1rpUMI=top4)
# 
# # generate random numbers according to the fitted distributions
# r.beta.ag <- c(rep(0, round(nsim*p0[1])), rbeta(nsim-round(nsim*p0[1]), shape1[1], shape2[1]))
# r.beta.ga <- c(rep(0, round(nsim*p0[2])), rbeta(nsim-round(nsim*p0[2]), shape1[2], shape2[2]))
# r.beta.ct <- c(rep(0, round(nsim*p0[3])), rbeta(nsim-round(nsim*p0[3]), shape1[3], shape2[3]))
# r.beta.tc <- c(rep(0, round(nsim*p0[4])), rbeta(nsim-round(nsim*p0[4]), shape1[4], shape2[4]))
# 
# # evaluate the fit -- QQplot
# png('QQplot.beta.include_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   r <- eval(parse(text=paste0('r.beta.', tolower(type)))) 
#   p <- eval(parse(text=paste0('p', type)))
#   qqplot(p, r, xlab='observed', ylab='fitted', main=type)
#   abline(a=0, b=1, col='red', lty=2)
# }
# dev.off()
# 
# # evaluate the fit -- density function in the tail
# png('tail.beta.include_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   p <- eval(parse(text=paste0('p', type)))
#   plot(density(p[p>0]), xlim=c(.0008, .002), ylim=c(0, 100), main=type, xlab='')
#   fit <- eval(parse(text=paste0('fit.beta.mle.', tolower(type))))
#   a <- fit$estimate[1]
#   b <- fit$estimate[2]
#   curve(dbeta(x, a, b), add=T, col='red', lty=2)
# }
# dev.off()
# 
# # evaluate the fit -- ecdf
# png('ECDF.beta.include_1rpUMI.png', height=1200, width=1200, res=200)
# par(mfrow=c(2,2))
# for(type in c('AG', 'GA', 'CT', 'TC')){
#   r <- eval(parse(text=paste0('r.beta.', tolower(type)))) 
#   p <- eval(parse(text=paste0('p', type)))
#   plot(ecdf(p), xlab='', main=type, ylim=c(0.95, 1))
#   lines(ecdf(r), lty=2, col='red')
#   legend(x='bottomright', legend=c('observed', 'fitted'), lty=1:2, col=c('black', 'red'))
# }
# dev.off()



#####################################
# paretoQPlot(pTC)
x0.AG <- 0.0025
x0.GA <- 0.0025
x0.CT <- 0.0020
x0.TC <- 0.0020

ul.AG <- 0.040
ul.GA <- 0.025
ul.CT <- 0.020
ul.TC <- 0.030

x0s <- thetas <- shape1s <- shape2s <- pBetas <- rep(NA, 4)

# fit beta distribution; 0's are imputed
for(i in 1:4){
  type <- types4[i]
  # fit pareto tails
  ul <- eval(parse(text=paste0('ul.', type)))
  p.tmp <- eval(parse(text=paste0('p', type)))
  x0 <- eval(parse(text=paste0('x0.', type)))
  pareto.theta <- paretoTail(p.tmp[p.tmp < ul], x0 = x0)$theta
  # fit beta distribution
  p.min <- min(p.tmp[p.tmp > 0])
  p.tmp1 <- ifelse(p.tmp == 0, runif(1, 0, p.min), p.tmp)
  fit <- fitdist(p.tmp1, distr='beta', method='mle')
  shape1 <- fit$estimate[1]
  shape2 <- fit$estimate[2]
  
  # simulate beta's
  tmp <- rbeta(5*nsim, shape1, shape2)
  tmp <- tmp[tmp <= x0]
  p.beta <- sum(p.tmp <= x0)/length(p.tmp)
  r.beta <- sample(tmp, round(nsim*p.beta))
  
  # simulate pareto's; remove extremely large data points
  r.pare <- x0 * runif(round(nsim * (1-p.beta)))^(-1/pareto.theta)
  p.max <- max(p.tmp)
  r.pare <- ifelse(r.pare <= ul, r.pare, runif(1, 0, p.max))
  assign(paste0('r.', tolower(type)), c(r.beta, r.pare)) 
  
  # save parameters
  x0s[i] <- x0
  thetas[i] <- pareto.theta
  shape1s[i] <- shape1
  shape2s[i] <- shape2
  pBetas[i] <- p.beta
}

parameters.include_1rpUMI <- data.frame(types4, x0s, thetas, shape1s, shape2s, pBetas)
colnames(parameters.include_1rpUMI) <- c('type', 'x0', 'theta', 'shape1', 'shape2', 'pBeta')
r.inc1.ag <- r.ag
r.inc1.ga <- r.ga
r.inc1.ct <- r.ct
r.inc1.tc <- r.tc

# evaluate the fit -- QQplot
png('QQplot.piecewise.include_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  r <- eval(parse(text=paste0('r.', tolower(type)))) 
  p <- eval(parse(text=paste0('p', type)))
  qqplot(p, r, xlab='observed', ylab='fitted', main=type)
  abline(a=0, b=1, col='red', lty=2)
}
dev.off()

# evaluate the fit -- full density
png('density.piecewise.include_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  p <- eval(parse(text=paste0('p', type)))
  plot(density(p), xlim=c(0, .002), ylim=c(0, 5000), main=type, xlab='')
  
  fit <- eval(parse(text=paste0('r.', tolower(type))))
  lines(density(fit), col='red', lty=2)
}
dev.off()

# evaluate the fit -- density function in the tail
png('tail.piecewise.include_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  p <- eval(parse(text=paste0('p', type)))
  plot(density(p), xlim=c(0.0008, .002), ylim=c(0, 100), main=type, xlab='')
  
  fit <- eval(parse(text=paste0('r.', tolower(type))))
  lines(density(fit), col='red', lty=2)
}
dev.off()

# evaluate the fit -- ecdf
png('ECDF.piecewise.include_1rpUMI.png', height=1200, width=1200, res=200)
par(mfrow=c(2,2))
for(type in c('AG', 'GA', 'CT', 'TC')){
  r <- eval(parse(text=paste0('r.', tolower(type)))) 
  p <- eval(parse(text=paste0('p', type)))
  plot(ecdf(p), xlab='', main=type, ylim=c(0.95, 1))
  lines(ecdf(r), lty=2, col='red')
  legend(x='bottomright', legend=c('observed', 'fitted'), lty=1:2, col=c('black', 'red'))
}
dev.off()

##################################################################
# save fit statistics for calculating p-values
##################################################################
bkg.error <- list(parameters.include_1rpUMI = parameters.include_1rpUMI, 
                  parameters.exclude_1rpUMI = parameters.exclude_1rpUMI, 
                  r.exc1.ag = r.exc1.ag, 
                  r.exc1.ga = r.exc1.ga,
                  r.exc1.ct = r.exc1.ct,
                  r.exc1.tc = r.exc1.tc,
                  r.inc1.ag = r.inc1.ag, 
                  r.inc1.ga = r.inc1.ga,
                  r.inc1.ct = r.inc1.ct,
                  r.inc1.tc = r.inc1.tc)

save(bkg.error, file='bkg.error.v2.4.RData')
#load(file='bkg.error.v2.4.RData')
