setwd('~/Desktop/svj_R_package stochastic intensity/')
options(digits=12)

set.seed(1)
seed = c(123, 321, 231, 132)

library(Rcpp)
library(BH)
library(nloptr)
library(ggplot2)
library(gridExtra)

source('conv_plot.R')
sourceCpp('Simulation.cpp')
sourceCpp('nll.cpp', verbose = F)

par = c(0.05, 4.5394, 0.0439, 0.3038, -0.6974, 1.1308, 3, -0.2151, 0.15*0.15)
y = SimSVJ(par, 2001, log(100), par[3], 1/252, seed[1], seed[2], seed[3], seed[4])
f <- function(p) {
    return(nll(p, y))
}

start.time = Sys.time()
a = f(par)
Sys.time() - start.time
print(a)

start.time = Sys.time()

lb = c(-0.5, 0.001, 1e-6, 0.001, -0.99, 0, 0, -1, 1e-6)
ub = c(0.5, 100, 1, 1, 0.99, 10, 100, 0, 1)
opts <- list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel=1.0e-4, "maxeval"=10000, print_level = 0)

res1 = nloptr(x0 = par, eval_f = f, opts = opts, lb = lb, ub = ub)
print(res1)
Sys.time() - start.time
check.converge(f, res1$solution)

p = converge.plot(res1$solution)
print(p)
ggsave('~/Desktop/svj_R_package stochastic intensity/distribution.eps', p, width = 10, height = 7)