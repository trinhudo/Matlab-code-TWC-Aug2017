% SUM OF SUM, DPS VG-AF, ASYMPTOTIC
function out = SOSforAsymptoticDPSVGAF(lSRk,lRkDm,snravg,snrth,eta,rho,SETfull)
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
            2*snrth*subsum*lRkDm/eta/rho*...
            (1-snrth*subsum/(1-rho)/snravg)*...
            log(sqrt(snrth*lRkDm*subsum/eta/rho/snravg));
        out = out + ((-1)^(ii+1))*A;
    end
end

