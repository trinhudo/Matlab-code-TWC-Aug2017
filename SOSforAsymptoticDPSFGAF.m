% SUM OF SUM, DPS FG-AF, ASYMPTOTIC
function out = SOSforAsymptoticDPSFGAF(lSRk,lSRk_kappa,lRkDm,snravg,snrth,eta,rho,SETfull)
out = 0;
for ii = 1:length(SETfull)
    elements = nchoosek(SETfull,ii); % n_1 < ... < n_t
    [row,column] = size(elements);
    for rr = 1:row
        subsum = 0;
        for cc = 1:column
            subsum = subsum + lSRk(elements(rr,cc));
        end
        %
        A = snrth*subsum/(1-rho) - ...
            (1-snrth*subsum/(1-rho)/snravg)*...
            (2*snrth*lRkDm*subsum/eta/rho/(1-rho)+...
            2*snrth*lSRk_kappa*lRkDm*subsum/eta/rho/(1-rho)/snravg)*...
            log(sqrt(snrth*(snravg+lSRk_kappa)*lRkDm*subsum/eta/rho/(1-rho)/(snravg^2)));
        %
        out = out + ((-1)^(ii+1))*A;
    end
end

