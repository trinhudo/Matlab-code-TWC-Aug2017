% DPS FG-AF, SIMULATION
function DPS_FGAF_sim = DPS_FGAF_simulation...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
count = 0;
%% CHANNELS
for mm = 1:M
    hSDm(:,mm) = random('Rayleigh',sqrt(1/lSDm(mm)/2),[1,Sim_times]);
end
for kk = 1:K
    hSRk(:,kk) = random('Rayleigh',sqrt(1/lSRk(kk)/2),[1,Sim_times]);
    for mm = 1:M
        hRkDm(:,kk,mm) = random('Rayleigh',sqrt(1/lRkDm(kk,mm)/2),[1,Sim_times]);
    end
    %
    kappak(kk) = snravg+lSRk(kk);
end
%% SNRs
snrSDm = snravg*(abs(hSDm).^2); 
% Received SNR at destinations
[snrSDb(:,1),D] = max(snrSDm,[],2); 
% Find the best destination
%
snrSRk = (1-rho)*snravg*(abs(hSRk).^2); 
% Received SNR at relays
[snrSRb(:,1),R] = max(snrSRk,[],2); 
% Find the best relay (NOTE: partial relay selection)
for yy = 1:Sim_times
    hSRb(yy,1) = hSRk(yy,R(yy));
    hRbDb(yy,1,1) = hRkDm(yy,R(yy),D(yy));
    kappab(yy,1) = kappak(R(yy));
end
snrRbDb = eta*rho*snravg*(abs(hSRb).^2).*(abs(hRbDb).^2);
% SNR of the second hop with best relay and destination
for jj = 1:Sim_times
    snr_relay_path(jj,:) = snrSRb(jj,:).*snrRbDb(jj,:)./...
        (kappab(jj).*(abs(hSRb(jj,:)).^2)+snrRbDb(jj,:));
    % e2e SNR of relaying channel
end
snr_e2e = max(snrSDb,snr_relay_path);
% e2e SNR of the system (selection combining)
%% COUNT EVENTS
for zz = 1:Sim_times
    if (snr_e2e(zz) < snrth)
        count = count + 1;
    end
end
DPS_FGAF_sim = count/Sim_times;
end