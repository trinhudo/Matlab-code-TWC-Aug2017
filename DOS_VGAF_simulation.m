% DOS VG-AF, SIMULATION
function DOS_VGAF_sim = DOS_VGAF_simulation...
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
end
%% SNRs
snrSDm = snravg*(abs(hSDm).^2); % received SNR at Dm
[snrSDb(:,1),D] = max(snrSDm,[],2); % find the best destination
%
snrSRk = (1-rho)*snravg*(abs(hSRk).^2); % received SNR at Rk
%
for yy = 1:Sim_times
    hRkDb(yy,:,1) = hRkDm(yy,:,D(yy));
end
snrRkDb = eta*rho*snravg*(abs(hSRk).^2).*(abs(hRkDb).^2);
% end-to-end SNR
snr_relay_path = snrSRk.*snrRkDb./(snrSRk+snrRkDb+1);
[snr_relay_path_best(:,1),R] = max(snr_relay_path,[],2);
snr_e2e = max(snrSDb,snr_relay_path_best);
%% COUNT OUTAGE EVENTS
for zz = 1:Sim_times
    if (snr_e2e(zz) < snrth)
        count = count + 1;
    end
end
DOS_VGAF_sim = count/Sim_times;
end
