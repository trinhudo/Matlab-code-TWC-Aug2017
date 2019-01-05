% DOS FG-AF, EXACT
function DOS_FGAF_exact_out = DOS_FGAF_exact...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon);
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
    for kk = 1:K
        kappak = snravg+lSRk(kk);
        prod = prod*(1 - ...
            exp(-mu*lSRk(kk))* ...
            sqrt(4*kappak*snrth*lSRk(kk)*lRkDm(kk,mm)/alpha/beta)*...
            besselk(1,sqrt(4*kappak*snrth*lSRk(kk)*lRkDm(kk,mm)/alpha/beta)));
    end
    temp1 = temp1 + Theta*prod;
end
DOS_FGAF_exact_out =  Omega*temp1;
end