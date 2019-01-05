% DPS DF, EXACT and APPROXIMATE
function [DPS_DF_exact,DPS_DF_approx] = DPS_DF_exact_approx...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon);
set_dest_full = 1:M;
set_relay_full = 1:K;
%
alpha = (1-rho)*snravg;
beta = eta*rho*snravg;
%% CALCULATIONS
Omega = 1 + SOSforOmega(lSDm,snrth,snravg,set_dest_full);
%
Theta = 0;
temp1 = 0;
temp3 = 0;
temp1_approx = 0;
temp3_approx = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp2 = 0;
    temp2_approx = 0;
    for kk = 1:K
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        temp3= SOSforExactDPSDF(lSRk,lRkDm(kk,mm),snrth,alpha,beta,set_relay_full);
        temp2 = temp2 + Xi*temp3;
        %
        temp3_approx= SOSforApproxDPSDF(lSRk,lRkDm(kk,mm),...
            snrth,alpha,beta,set_relay_full);
        temp2_approx = temp2_approx + Xi*temp3_approx;
        %
    end
    temp1 = temp1 + Theta*temp2;
    %
    temp1_approx = temp1_approx + Theta*temp2_approx;
end
%
DPS_DF_exact = Omega*(1-temp1);
DPS_DF_approx = Omega*(1-temp1_approx);
end