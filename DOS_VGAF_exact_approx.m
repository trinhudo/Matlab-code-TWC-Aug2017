% DOS VG-AF, EXACT and APPROXIMATE
function [DOS_VGAF_exact,DOS_VGAF_approx] = DOS_VGAF_exact_approx...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%
alpha = (1-rho)*snravg;
beta = eta*rho*snravg;
mu = snrth/alpha;
%% CALCULATIONS
Omega = 1 + SOSforOmega(lSDm,snrth,snravg,set_dest_full);
%
temp1 = 0;
temp1_approx = 0;
Theta = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    prod = 1;
    prod_approx = 1;
    for kk = 1:K
        Y = @(y) exp(- y*lSRk(kk) - ...
            (alpha.*snrth.*y+snrth).*lRkDm(kk,mm)./...
            ((alpha.*y-snrth).*beta.*y));
        temp2 = integral(Y,mu,inf);
        prod = prod*(1 - lSRk(kk)*temp2);
        %
        prod_approx = prod_approx*(1 - exp(-mu*lSRk(kk))* ...
            sqrt(4*snrth*lSRk(kk)*lRkDm(kk,mm)/beta)*...
            besselk(1,sqrt(4*snrth*lSRk(kk)*lRkDm(kk,mm)/beta)));
    end
    temp1 = temp1 + Theta*prod;
    temp1_approx = temp1_approx + Theta*prod_approx;
end
DOS_VGAF_exact =  Omega*temp1;
DOS_VGAF_approx =  Omega*temp1_approx;
end