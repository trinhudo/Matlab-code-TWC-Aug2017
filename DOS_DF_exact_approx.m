% DOS DF, EXACT AND APPROXIMATE
function [DOS_DF_exact,DOS_DF_approx] = DOS_DF_exact_approx...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%
alpha = (1-rho)*snravg;
beta = eta*rho*snravg;
mu = snrth/alpha;
%
Omega = 1 + SOSforOmega(lSDm,snrth,snravg,set_dest_full);
%
temp1 = 0;
temp1_approx = 0;
Theta = 0;
%% CALCULATIONS
for mm = 1:M
    set_dest_minus = setdiff(set_dest_full,mm);
    Theta = (1 + SOSforThetaXi(lSDm(mm),lSDm,set_dest_minus));
    prod = 1;
    prod_approx = 1;
    for kk = 1:K
        Y = @(y) exp(- y*lSRk(kk) - snrth./beta./y.*lRkDm(kk,mm));
        temp2 = integral(Y,mu,inf);
        prod = prod*(1 - lSRk(kk)*temp2);
        %
        % New approx. (2016-08-20)
        nu = lSRk(kk);
        xi = snrth*lRkDm(kk,mm)/beta;
        %
        a = mu;
        b = xi;
        c = nu;
        %
        A1 = exp(-a*c)/c;
        A2 = -b*igamma(0,a*c);
        A3 = 0;
        for nn=2:6
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
        prod_approx = prod_approx*(1 - lSRk(kk)*NewApprox);
    end
    temp1 = temp1 + Theta*prod;
    temp1_approx = temp1_approx + Theta*prod_approx;
end
DOS_DF_exact =  Omega*temp1;
DOS_DF_approx =  Omega*temp1_approx;
end
