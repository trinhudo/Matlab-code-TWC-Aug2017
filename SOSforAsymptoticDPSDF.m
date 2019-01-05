% SUM OF SUM, DPS DF, ASYMPTOTIC
function out = SOSforAsymptoticDPSDF(lSRk,lRkDm,snravg,snrth,eta,rho,SETfull)
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
        A = snrth*subsum/(1-rho) -...
            snrth*lRkDm*subsum/eta/rho*...
            (0.577215 + log(snrth*subsum/(1-rho)/snravg));
        %
        out = out + ((-1)^(ii+1))*A;
    end
end

