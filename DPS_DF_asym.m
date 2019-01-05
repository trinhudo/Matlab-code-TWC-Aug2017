% DPS DF, ASYMPTOTIC
function DPS_DF_asym_out = DPS_DF_asym...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%% CALCULATIONS
prod_omega = 1;
for pp = 1:M
    prod_omega = prod_omega*snrth*lSDm(pp);
end
%
Omega_asym = ((1/snravg)^(M+1))*prod_omega;
%
Theta = 0;
temp1_asym = 0;
temp3_asym = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp2_asym = 0;
    for kk = 1:K
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        %
        temp3_asym= ...
            SOSforAsymptoticDPSDF(lSRk,lRkDm(kk,mm),snravg,snrth,eta,rho,set_relay_full);
        temp2_asym = temp2_asym + Xi*temp3_asym;
        %
    end
    %
    temp1_asym = temp1_asym + Theta*temp2_asym;
end
%
DPS_DF_asym_out = Omega_asym*temp1_asym;
end