% DPS VG-AF, EXACT
function DPS_FGAF_exact_out = DPS_FGAF_exact...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%
alpha = (1-rho)*snravg;
beta = eta*rho*snravg;
%% CALCULATIONS
Omega = 1 + SOSforOmega(lSDm,snrth,snravg,set_dest_full);
%
temp3 = 0;
Theta_m = 0;
temp1 = 0;
A2approx = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp2 = 0;
    for kk = 1:K
        kappak = snravg+lSRk(kk);
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        temp1 = 1 - ...
            SOSforExactDPSFGAF(lSRk,lRkDm(kk,mm),snrth,...
            alpha,beta,kappak,set_relay_full);
        temp2 = temp2 + Xi*temp1;
        %
    end
    temp3 = temp3 + Theta*temp2;
end
%
DPS_FGAF_exact_out = Omega*temp3;
end
