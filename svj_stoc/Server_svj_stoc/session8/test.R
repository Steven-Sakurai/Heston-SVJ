set.seed(515)

n = 50
len = 1001
require(Rcpp)
require(BH)
require(nloptr)
require(ggplot2)
require(gridExtra)
options(digits=8)
source('../conv_plot.R')
sourceCpp('../Simulation.cpp')
par = c(0.05, 3.9, 0.08, 0.3038, -0.6974, 3.2, -0.3551, 0.0967*0.0967)
sourceCpp('../nll.cpp', verbose = F)

opts <- list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel = 1.0e-6, maxeval = 10000)
# set lower and upper bound
lb = c(-0.5, 0.001, 1e-6, 0.001, -0.99, 0, -1, 1e-8)
ub = c(0.5, 100, 1, 1, 0.99, 10, 0, 1)

# run n paths
converge.list = vector(mode="numeric", length = n)
est.par.list = data.frame(mu = double(), kappa = double(), theta = double(), xi = double(), rho = double(), lambda = double(), mu_s = double(), sigma_s = double())
i = 1
while(i <= n) {
    seed.list = sample(10000, 4)
    y = SimSVJ(par, len, log(100), par[3], 1/252, seed.list[1], seed.list[2], seed.list[3], seed.list[4])
	f <- function(p) {
	    return(nll(p, y))
	}
    start = Sys.time()
    result = nloptr(x0 = par, eval_f = f, 
       opts = opts, 
       lb = lb, ub = ub)
    if(result$status < 0)
        next;
	est.res = as.numeric(result$solution)
	conv.res = check.converge(f, result$solution)
	write.table(est.res, file = "./est.csv", sep = ",", append = T, row.names = F, col.names = F, quote = F)
	write.table(conv.res, file = "./conv.csv", sep = ",", append = T, row.names = F, col.names = F, quote = F)
    print(i)
    print(Sys.time() - start)
    print(est.res)
    est.par.list[i, ] = as.numeric(est.res)
    converge.list[i] = conv.res
	i = i + 1
}

write.csv(file = "convResult.csv", converge.list)
write.csv(file = "estResult.csv", est.par.list)
