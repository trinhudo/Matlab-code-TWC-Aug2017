% DPS FG-AF, ASYMPTOTIC
function DPS_FGAF_asym_out = DPS_FGAF_asym...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%% CALCULATIONS
% Omega
prod_omega = 1;
for pp = 1:M
    prod_omega = prod_omega*snrth*lSDm(pp);
end
%
Omega_asym = ((1/snravg)^(M+1))*prod_omega;%
%
temp3_asym = 0;
temp1_asym = 0;
Theta = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp2_asym = 0;
    for kk = 1:K
        lSRk_kappa = lSRk(kk);
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        temp1_asym = SOSforAsymptoticDPSFGAF(lSRk,lSRk_kappa,lRkDm(kk,mm),snravg,snrth,eta,rho,set_relay_full);
        temp2_asym = temp2_asym + Xi*temp1_asym;
        %
    end
    temp3_asym = temp3_asym + Theta*temp2_asym;
end
%
DPS_FGAF_asym_out = Omega_asym*temp3_asym;
end
