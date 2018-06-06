clc
clear all
close all

%%% LOAD DURING INJECTION
load lamda.mat;

load S_LP_DI;
S_UPD_DI=S_LP_DI.Supd;



t0=6; % post injection starting time
tend=15; % end time for calculation
dt=0.25; % time interval 6 hours = 0.25 days
tp=t0:dt:tend; % post injection time vector 

% observed rates during 15 days
load rate_obs15d.mat;
obs=rate_obs15d(1:60);
% observed rates during injection (di)
obs_di=obs(1:23); 
% observed rates post injection (pi)
obs_pi=obs(24:59);
% small change for numerical comp. reasons
for ig=1:length(obs_pi);
    if obs_pi(ig)==0;
        obs_pi(ig)=0.1;
    end
end


t00=5.75; % injection stopping time

%% Estimation of p, updating every 6 hours
xx=log10((mean(obs_di)./obs_pi)); xx=xx';
yy=log10(tp(2:end)./t00);
p_up=xx./yy;
p_upd(1)=2; % initial p value
p_upd(2:length(p_up)+1)=p_up;


% POST-INJECTION UPDATED AS MORE DATA ARRIVES
S_PI_upd=mean(obs_di)./((tp./t00).^p_upd);

dt=0.25;
tmax=1.0:dt:5.5;
pp=4;
for i=1:length(tmax);
    S(i).pi=mean(obs_di(1:(tmax(i)/dt)))./((tp./t00).^pp);
end


time=0.25:0.25:15;
%% pi + di
Supd=[S_UPD_DI' S_PI_upd];
S1D=[lamda(1).lamda_sr S(1).pi];
S2D=[lamda(5).lamda_sr S(5).pi];
S3D=[lamda(9).lamda_sr S(9).pi];
S4D=[lamda(13).lamda_sr S(13).pi];
S5D=[lamda(17).lamda_sr S(17).pi];


figure;
plot(time,obs,'r','LineWidth',8); hold on;
plot(time,Supd,'Color',[0.7 0.7 0.7],'LineWidth',10); hold on;
plot(time,S1D,'-','Color',[0.4 0.4 0.4],'LineWidth',7); hold on;
plot(time,S2D,'-','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(time,S3D,'b','LineWidth',2); hold on;
plot(time,S4D,'c','LineWidth',2); hold on;
% plot(time,S5D,'g','LineWidth',2); hold on;
legend('obs','SR','1 day','2 days','3 days','4 days');
ylabel('Seismicity rates','FontSize',14);
xlabel('Time (days)','FontSize',14);
grid on;





lp=1:dt:4.0;
for j=1:length(lp);
    SLP = [lamda(j).lamda_sr S(j).pi];
    SR_LP(j).MODELS=SLP;
    for g=1:length(obs);
        log_like(g) = sum(log(poisspdf(obs(g),SLP(g))));
        log_like_upd(g) = sum(log(poisspdf(obs(g),Supd(g))));
    end
    lls(j)=sum(log_like(17:end));
    
    clear loglike
end

lls_upd=sum(log_like_upd(17:end));
lls_up(1:length(lp))=lls_upd;




load llse.mat;
load lls_SR.mat;

t=1:0.25:4.0;
figure;

plot(t,llse(1:length(t)),'bo-','LineWidth',2); hold on;
plot(lp,lls,'ko-','LineWidth',2); hold on;
plot(lp,lls_up,'k--','LineWidth',6); hold on;
plot(t,lls_SR(1:length(t)),'b--','LineWidth',6);
legend('E5','SR','E5-reference','SR-reference');
xlabel('Learning Period (days)','FontSize',18);
ylabel('log-likelihood','FontSize',18);
grid on;



loglike_LP.E5=llse(1:length(t));
loglike_LP.SR=lls;
loglike_LP.E5ref=lls_up;
loglike_LP.SRref=lls_SR(1:length(t));

% savefile='loglike_LP.mat';
% save(savefile,'loglike_LP');



% figure;
% for ii=1:13;
%     plot(SR_LP(ii).MODELS); hold on;
% end
% savefile='SR_LP.mat';
% save(savefile,'SR_LP');








