% OMEGA
% sum-of-sum function for m = 1 to M (full case)
% K = 5; % number of relays
% M = 5; % number of destinations
% beta = 2.7; % path-loss exponent
% % lambda (average channel gains)
% lSD = lambda(K,M,beta);
% Rth = 2; % bits/s/Hz
% snrth = 2.^(2.*Rth) - 1;
% snravg_dB = 10; % P/N0 in dB
% snravg = 10.^(snravg_dB./10);
% num_dest = M;
% SET = 1:num_dest;
%
function out = SOSforOmega(lSD,snrth,snravg,SETfull)
out = 0;
for ii = 1:length(SETfull)
    elements = nchoosek(SETfull,ii); % n_1 < ... < n_t
    [row,column] = size(elements);
    for rr = 1:row
        subsum = 0; 
        for cc = 1:column
            subsum = subsum + lSD(elements(rr,cc));
        end
        out = out + ((-1)^(ii)) * exp(-snrth/snravg*subsum);
    end
end