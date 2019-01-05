% DOS and DPS, 3 SCHEMES, EXACT and ASYMPTOTIC
close all
clear all
%% PARAMETERS
K = 3; % number of relays
M = 3; % number of destinations
rho = 0.5; % power splitting ratio
snravg_dB = -10:1:20; % in dB
snravg = 10.^(snravg_dB./10);
Rth = 1; % bits/s/Hz
snrth = 2.^(2.*Rth) - 1;
espsilon = 2.7; % path loss exponent
eta = 0.7; % energy conversion efficiency
%
[lSDm,lSRk,lRkDm] = lambda(K,M,espsilon); % lambda
set_dest_full = 1:M;
set_relay_full = 1:K;
%
Sim_times = 10^0; % for Monte-Carlo
%% SIMULATIONS
for ss = 1:length(snravg_dB)
    disp(['SNR=' num2str(snravg_dB(ss)) 'dB']);
    % DOS DF
    [DOS_DF_exact_out(ss),DOS_DF_approx_out(ss)] = DOS_DF_exact_approx...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DOS_DF_asym_out(ss) = DOS_DF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    % DOS VG-AF
    [DOS_VGAF_exact_out(ss),DOS_VGAF_approx_out(ss)] = DOS_VGAF_exact_approx...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DOS_VGAF_asym_out(ss) = DOS_VGAF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    % DOS FG-AF
    DOS_FGAF_exact_out(ss) = DOS_FGAF_exact...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DOS_FGAF_asym_out(ss) = DOS_FGAF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    %%
    % DPS DF
    [DPS_DF_exact_out(ss),DPS_DF_approx(ss)] = ...
        DPS_DF_exact_approx...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DPS_DF_asym_out(ss) = DPS_DF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    % DPS VG-AF
    [DPS_VGAF_exact_out(ss),DPS_VGAF_approx(ss)] = DPS_VGAF_exact_approx...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DPS_VGAF_asym_out(ss) = DPS_VGAF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    % DPS FG-AF
    DPS_FGAF_exact_out(ss) = DPS_FGAF_exact...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
    DPS_FGAF_asym_out(ss) = DPS_FGAF_asym...
        (K,M,rho,snrth,snravg(ss),espsilon,eta,Sim_times);
end

%% PLOT
semilogy(snravg_dB,DOS_DF_exact_out,'bs-',...
    'MarkerFaceColor','b')
hold on
semilogy(snravg_dB,DOS_VGAF_exact_out,'b^-',...
    'MarkerFaceColor','b')
hold on
semilogy(snravg_dB,DOS_FGAF_exact_out,'bo-',...
    'MarkerFaceColor','b')
hold on
semilogy(snravg_dB,DPS_DF_exact_out,'rs-')
hold on
semilogy(snravg_dB,DPS_VGAF_exact_out,'r^-')
hold on
semilogy(snravg_dB,DPS_FGAF_exact_out,'ro-')
hold on
%
semilogy(snravg_dB,DOS_DF_asym_out,'b:')
hold on
semilogy(snravg_dB,DOS_VGAF_asym_out,'b:')
hold on
semilogy(snravg_dB,DOS_FGAF_asym_out,'b:')
hold on
semilogy(snravg_dB,DPS_DF_asym_out,'r:')
hold on
semilogy(snravg_dB,DPS_VGAF_asym_out,'r:')
hold on
semilogy(snravg_dB,DPS_FGAF_asym_out,'r:')
%%
xlabel('Transmit SNR (dBm)')
ylabel('Outage Probability')
legend('DOS DF','DOS VG-AF','DOS FG-AF',...
    'DPS DF','DPS VG-AF','DPS FG-AF',...
    'Location','SouthWest')
%
set(gca,'XTick',-10:10:20)
set(gca,'XTickLabel',[20 30 40 50] );
axis([-10 20 10^-8 1])