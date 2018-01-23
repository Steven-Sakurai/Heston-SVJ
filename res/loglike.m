% Q: why calculate the full filter moment? 




function [MDensity, logPy, FMoment] = loglike(Par, YObs, N, n)

tic

mu = Par(1);
kappa = Par(2);
theta = Par(3);
xi = Par(4);
rho = Par(5);

nbasis = 21;
Basis = zeros(N + 1, 4*nbasis);
FMoment = zeros(N + 1, 4);
Py = zeros(N, 1);
logPy = zeros(N, 1);
y = YObs;

Basis(1, :) = InitialMoments(Par); % M_k
%filter moment, sum over truncated filter moment, 0th~3rd order
FMoment(1, 1) = sum(Basis(1, 1 : nbasis)); 
FMoment(1, 2) = sum(Basis(1, nbasis + 1 : 2*nbasis));
FMoment(1, 3) = sum(Basis(1, 2*nbasis + 1 : 3*nbasis));
FMoment(1, 4) = sum(Basis(1, 3*nbasis + 1 : 4*nbasis));

for j = 1 : N %online induction
    
	% 重点
	[MarTranCoef, BasisCoef] = BasisCoefficients_cpp_EA_Heston(Par, y(j), y(j + 1));
	% 重点
	
	Py(j) = MarTranCoef*Basis(j, :)'; 
	% L_i = \sum(k = 1-84) \alpha_k M_k
    logPy(j)=log(Py(j));
    
	Basis(j + 1, :) = transpose(BasisCoef*Basis(j, :)'/Py(j));
	
    FMoment(j + 1, 1) = sum(Basis(j + 1, 1 : nbasis));
    FMoment(j + 1, 2) = sum(Basis(j + 1, nbasis + 1 : 2*nbasis));
    FMoment(j + 1, 3) = sum(Basis(j + 1, 2*nbasis + 1 : 3*nbasis));
    FMoment(j + 1, 4) = sum(Basis(j + 1, 3*nbasis + 1 : 4*nbasis));
end

MDensity = -sum(logPy); 

if ~(isreal(MDensity))
	MDensity = 0;
end


fprintf('%d %f %f %f %f %f %f\n',n,mu,kappa,theta,xi,rho,MDensity)

end
