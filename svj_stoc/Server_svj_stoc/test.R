
n = 50
len = 1001
par = c(0.05, 2, 0.08, 0.5038, -0.7, 5, 50, -0.3551, 0.15)

require(Rcpp)
require(BH)
require(nloptr)
require(ggplot2)
require(gridExtra)
options(digits=8)
source('../conv_plot.R')
sourceCpp('../Simulation.cpp')
sourceCpp('../nll.cpp', verbose = F)

opts <- list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel = 1.0e-4, maxeval = 10000)
# set lower and upper bound
lb = c(-0.5, 0.001, 1e-6, 0.001, -0.99, 0, 0, -1, 0)
ub = c(1, 100, 1, 1, 0.99, 100, 1000, 0, 1)

npar = 9
# run n paths
converge.list = matrix(, nrow = n, ncol = npar)
est.par.list = matrix(, nrow = n, ncol = npar)
success.list = vector(mode="numeric", length = n)

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
	est_list = as.numeric(result$solution)
	c_list = check.converge(f, result$solution)
	m_list = check.mono(f, result$solution)
    flag = all(m_list | c_list)

    write.table(est_list, file = "./est.csv", sep = ",", append = T, row.names = F, col.names = F, quote = F)
	write.table(c_list, file = "./conv.csv", sep = ",", append = T, row.names = F, col.names = F, quote = F)
    write.table(flag, file = "./success.csv", sep = ",", append = T, row.names = F, col.names = F, quote = F)
    
    print(i)
    print(Sys.time() - start)
    print(est_list)
    
    est.par.list[i, ] = as.numeric(est_list)
    converge.list[i, ] = c_list
    success.list[i] = flag
	i = i + 1
}

write.csv(file = "convResult.csv", converge.list)
write.csv(file = "estResult.csv", est.par.list)
write.csv(file = "successResult.csv", success.list)
