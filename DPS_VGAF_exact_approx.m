% DPS VG-AF, EXACT and SIMULATION
function [DPS_VGAF_exact,DPS_VGAF_approx] = DPS_VGAF_exact_approx...
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
temp3 = 0;
temp_approx_3 = 0;
Theta = 0;
temp1 = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp2 = 0;
    temp_approx_2 = 0;
    for kk = 1:K
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        temp1 = 1 - SOSforExactDPSVGAF(lSRk,lRkDm(kk,mm),snrth,...
            alpha,beta,set_relay_full);
        temp2 = temp2 + Xi*temp1;
        %
        temp_approx_1 = 1 - SOSforApproxDPSVGAF...
            (lSRk,lRkDm(kk,mm),snrth,alpha,beta,set_relay_full);
        temp_approx_2 = temp_approx_2 + Xi*temp_approx_1;
    end
    temp3 = temp3 + Theta*temp2;
    temp_approx_3 = temp_approx_3 + Theta*temp_approx_2;
end
%
DPS_VGAF_exact = Omega*temp3;
DPS_VGAF_approx = Omega*temp_approx_3;
end
