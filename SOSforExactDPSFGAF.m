% SUM OF SUM, DPS FG-AF, EXACT
function out = SOSforExactDPSFGAF(lSRk,lRkDm,snrth,alpha,beta,kappa,SETfull)
out = 0;
for ii = 1:length(SETfull)
    elements = nchoosek(SETfull,ii); % n_1 < ... < n_t
    [row,column] = size(elements);
    for rr = 1:row
        subsum = 0;
        for cc = 1:column
            subsum = subsum + lSRk(elements(rr,cc));
        end
        out = out + ((-1)^(ii+1)) * exp(-snrth*subsum/alpha)*...
            sqrt(4*kappa*snrth*lRkDm*subsum/alpha/beta)*...
            besselk(1,sqrt(4*kappa*snrth*lRkDm*subsum/alpha/beta));
    end
end