% DOS VG-AF, ASYMPTOTIC
function DOS_VGAF_asym_out = DOS_VGAF_asym...
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
Omega_asym = ((1/snravg)^(M+K))*prod_omega;
%
temp1_asym = 0;
Theta = 0;
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    prod_asym = 1;
    for kk = 1:K
        %
        prod_asym = prod_asym*(snrth*lSRk(kk)/(1-rho) - ...
            2*snrth*lSRk(kk)*lRkDm(kk,mm)/eta/rho*...
            (1-snrth*lSRk(kk)/(1-rho)/snravg)*...
            log(sqrt(snrth*lSRk(kk)*lRkDm(kk,mm)/eta/rho/snravg)));
    end
    temp1_asym = temp1_asym + Theta*prod_asym;
end
DOS_VGAF_asym_out =  Omega_asym*temp1_asym;
end