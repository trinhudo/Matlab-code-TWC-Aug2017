% SUM OF SUM, DPS DF, APPROXIMATION
function out = SOSforApproxDPSDF(lSRk,lRkDm,snrth,alpha,beta,SETfull)
out = 0;
mu = snrth/alpha;
for ii = 1:length(SETfull)
    elements = nchoosek(SETfull,ii); % n_1 < ... < n_t
    [row,column] = size(elements);
    for rr = 1:row
        subsum = 0;
        for cc = 1:column
            subsum = subsum + lSRk(elements(rr,cc));
        end 
        % New approx. (2016-08-20)
        nu = subsum;
        xi = snrth*lRkDm/beta;
        % 
        a = mu;
        b = xi;
        c = nu;
        %
        A1 = exp(-a*c)/c;
        A2 = -b*igamma(0,a*c);
        A3 = 0;
        for nn=2:9
            B1 = ((-1)^nn)*(b^nn)/(factorial(nn));
            B21 = exp(-(a*c));
            B22 = 0;
            for vv = 1:(nn-1)
                temp = (factorial(vv-1))*((-c)^(nn-vv-1))/...
                    ((factorial(nn-1))*(a^vv));
                B22 = B22 + temp;
            end
            B23 = ((-c)^(nn-1))/(factorial(nn-1))*(ei(-a*c));
            B2 = B21*B22-B23;
            A3 = A3 + B1*B2;
        end
        NewApprox = A1+A2+A3;
        %
        out = out + ((-1)^(ii+1)) * subsum * NewApprox;
%         out = out + ((-1)^(ii+1)) * subsum * ...
%             ((exp(-mu*subsum))/subsum -...
%             snrth/beta*lRkDm*igamma(0,mu*subsum));
    end
end

