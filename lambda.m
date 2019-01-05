% calculating lambda based on distance and path-loss
function [lSDm,lSRk,lRkDm] = lambda(K,M,beta)
dSD = 10;
dSR = 3;
dRD = 7;

lSDm = (1e-3)*(dSD^beta)*ones(1,M);
lSRk = (1e-3)*(dSR^beta)*ones(1,K);
lRkDm = (1e-3)*(dRD^beta)*ones(K,M);

% % clear all
% % close all
% % K = 5;
% % M = 5;
% % beta = 2.7;
% S   = [0, .5];
% % Asymmetric topology
% D1  = [1,.45];
% D2  = [1,.5];
% D3  = [1,.55];
% D4  = [1,.4];
% D5  = [1,.6];
% R1  = [.4,.45];
% R2  = [.4,.5];
% R3  = [.4,.55];
% R4  = [.4,.4];
% R5  = [.4,.6];
% %
% dest = [D1;D2;D3;D4;D5];
% relay = [R1;R2;R3;R4;R5];
% % Calculate average channel gain \lambda
% lSDm = zeros(1,M);
% lSRk = zeros(1,K);
% lRkDm = zeros(K,M);
% for mm = 1:M
%     lSDm(mm) = (sqrt(abs(S(1) - dest(mm,1))^2 + ...
%         abs(S(2)-dest(mm,2))^2))^beta;
% end
% for kk = 1:K
%     lSRk(kk) = (sqrt(abs(S(1) - relay(kk,1))^2 + ...
%         abs(S(2)-dest(kk,2))^2))^beta;
%     for mm = 1:M
%         lRkDm(kk,mm) = (sqrt(abs(relay(kk,1) - dest(mm,1))^2 + ...
%             abs(relay(kk,2) - dest(mm,2))^2))^beta;
%     end
% end
% % %% Plot
% % figure(1)
% % % Topology
% % plot(S(1),S(2),'ks','MarkerSize',10,'MarkerFaceColor','b');
% % hold on
% % for mm = 1:M
% %     plot(dest(mm,1),dest(mm,2),'ks','MarkerSize',10,'MarkerFaceColor','r');
% %     hold on
% % end
% % for kk = 1:K
% %     plot(relay(kk,1),relay(kk,2),'ko','MarkerSize',10,'MarkerFaceColor','k');
% %     hold on
% %     plot([S(1),relay(kk,1)],[S(2),relay(kk,2)],'b-');
% %     plot([S(1),dest(kk,1)],[S(2),dest(kk,2)],'r-');
% %     hold on
% %     for mm = 1:M
% %         plot([relay(kk,1),dest(mm,1)],[relay(kk,2),dest(mm,2)],'--');
% %         hold on
% %     end
% % end
% % %
% % axis([0 1 0 1])
% % pbaspect([1 1 1])