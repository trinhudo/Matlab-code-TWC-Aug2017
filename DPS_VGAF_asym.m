% DPS VG-AF, ASYMPTOTIC
function DPS_VGAF_asym_out = DPS_VGAF_asym...
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
Omega_asym = ((1/snravg)^(M+1))*prod_omega;
%
temp_asym_3 = 0;
Theta = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    Xi = 0;
    temp_asym_2 = 0;
    for kk = 1:K
        set_relay_minus = setdiff(set_relay_full,kk);
        Xi = (1 + SOSforThetaXi(lSRk(kk),lSRk,set_relay_minus));
        temp_asym_1 = SOSforAsymptoticDPSVGAF(lSRk,lRkDm(kk,mm),snravg,snrth,eta,rho,set_relay_full);
        temp_asym_2 = temp_asym_2 + Xi*temp_asym_1;
    end
    temp_asym_3 = temp_asym_3 + Theta*temp_asym_2;
end
%
DPS_VGAF_asym_out = Omega_asym*temp_asym_3;
end
