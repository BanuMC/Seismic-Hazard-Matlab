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
% ylabel('Seismicity rates','FontSize',14);
% xlabel('Time (days)','FontSize',14);
grid on;
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',18,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within each time bin','FontSize',18,'FontName','Times')



Mmin=0.9;
Mmax=3.0;
dM=0.1;

%% Calculate rates for magnitude bins
[D15_forecast_magbin,observed_magbin,mag,D1_difr,D1_mb]=calc_rates_magbins2(Supd,Mmin,Mmax,dM);
[D1_forecast_magbin,observed_magbin,mag,D2_difr,D2_mb]=calc_rates_magbins2(S1D,Mmin,Mmax,dM);
[D2_forecast_magbin,observed_magbin,mag,D3_difr,D3_mb]=calc_rates_magbins2(S2D,Mmin,Mmax,dM);
[D3_forecast_magbin,observed_magbin,mag,D4_difr,D4_mb]=calc_rates_magbins2(S3D,Mmin,Mmax,dM);
[D4_forecast_magbin,observed_magbin,mag,SR_difr,SR_mb]=calc_rates_magbins2(S4D,Mmin,Mmax,dM);

mb2_SR=sum(SR_forecast_magbin(12:15))

% Compare forecasted and observed seismicity for each magnitude bin
figure;
plot(mag(2:end),observed_magbin,'r','LineWidth',8); hold on;
plot(mag(2:end),D15_forecast_magbin,'Color',[0.7 0.7 0.7],'LineWidth',10); hold on;grid on;
plot(mag(2:end),D1_forecast_magbin,'-','Color',[0.4 0.4 0.4],'LineWidth',7); hold on; 
plot(mag(2:end),D2_forecast_magbin,'-','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(mag(2:end),D3_forecast_magbin,'b','LineWidth',2); hold on;
plot(mag(2:end),D4_forecast_magbin,'c','LineWidth',2); hold on;
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M','FontSize',14); 
% legend('observed','SR','day1','day2','day3','day4');
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','M_w','FontSize',18,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within each magnitude bin','FontSize',18,'FontName','Times')







