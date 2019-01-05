% DOS FG-AF, SIMULATION
function DOS_FGAF_sim = DOS_FGAF_simulation...
    (K,M,rho,snrth,snravg,espsilon,eta,Sim_times)
%% PARAMETERS
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon);
set_dest_full = 1:M;
set_relay_full = 1:K;
%
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
snrSRk = (1-rho)*snravg*(abs(hSRk).^2);
% Received SNR at relays
for yy = 1:Sim_times
    hRkDb(yy,:,1) = hRkDm(yy,:,D(yy));
end
snrRkDb = eta*rho*snravg*(abs(hSRk).^2).*(abs(hRkDb).^2);
% SNR of the second hop with the given best destination
% end-to-end SNR
for jj = 1:Sim_times
    snr_relay_path(jj,:) = snrSRk(jj,:).*snrRkDb(jj,:)./...
        (kappak.*(abs(hSRk(jj,:)).^2) + snrRkDb(jj,:));
    % e2e SNR of the relaying channel
end
[snr_relay_path_best(:,1),R] = max(snr_relay_path,[],2);
% Find the best relay
snr_e2e = max(snrSDb,snr_relay_path_best);
% e2e SNR of the system
%% COUNT EVENTS
for zz = 1:Sim_times
    if (snr_e2e(zz) < snrth)
        count = count + 1;
    end
end
DOS_FGAF_sim = count/Sim_times;
end
