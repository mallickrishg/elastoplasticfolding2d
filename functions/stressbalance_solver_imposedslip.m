function slip = stressbalance_solver_imposedslip(taudinit,dtaud,taun,rcv,K,scf)
% stress balance solver
% solves the elasto-plastic governing equations for boundary elements
% |(Kd + mu*Kn)*ds| <= frictional resistance using a linear-inequality constraint solver
% taudinit - initial stress
% dtaud - incremental stress input into the system
% taun - fault normal stress (controld frictional strength)
% rcv - fault object (unicycle) - provide friction coefficients
% K - fault shear traction kernel
% scf - structure that contains the scaling factor 'cf' (convert from Pa -> MPa) to help solver converge
%       may also contain Ksigma - fault-normal traction kernel
% Rishav Mallick, EOS, 2020
% 
% converp to MPa to help calculations
if isstruct(scf)==0
	cf=scf;
else
	cf=scf.cf;
	Ksigma = scf.Ksigma./cf;
	option = scf.option;%1 - only tau; 2 - tau and sigma
end

options = optimoptions('lsqlin','MaxIterations',1e3);

tau = (taudinit + dtaud)./cf;
taun = taun./cf;
K = K./cf;

    
index = find(rcv.mu0>1);
Li = length(index);
index = rcv.mu0>1;

K = K(~index,~index);
tau = tau(~index);
taun = taun(~index);
mu0 = rcv.mu0(~index);

tausign = sign(tau);

N = rcv.N - Li;

A = [repmat(tausign,1,N).*K;...
    -repmat(tausign,1,N).*K];
b = [mu0.*taun - tausign.*tau;...
    tausign.*tau];
G = eye(N);
d = zeros(N,1);

Alin = [];
blin = [];

dslip = lsqlin(G,d,A,b,Alin,blin,[],[],[],options);

dslipvec = zeros(rcv.N,1);
dslipvec(~index) = dslip;

slip = dslipvec;





end