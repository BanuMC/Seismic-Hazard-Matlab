clc
clear all
close all

load COMB_MODEL_LP.mat;
load E5_LP.mat;
load SR_LP.mat;
load rate_obs15d.mat;

obs=rate_obs15d(1:60);
SR_D3=SR_LP(9).MODELS;
E5_D3=E5_LP(9).MODELS;
COM_D3=COMB_MODEL_LP.com9;


time=0.25:0.25:15;

figure;
plot(time,obs,'r','LineWidth',6); hold on;
plot(time,SR_D3,'-','Color',[0.4 0.4 0.4],'LineWidth',7); hold on;
plot(time,E5_D3,'-','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(time,COM_D3,'b','LineWidth',8); hold on;
legend('obs','SR','E5','COMB');
% ylabel('Seismicity rates','FontSize',14);
% xlabel('Time [days]','FontSize',14);
grid on;
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',18,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within six hours','FontSize',18,'FontName','Times')

Mmin=0.9;
Mmax=3.0;
dM=0.1;


%% Calculate rates for magnitude bins
[OBS_magbin,observed_magbin,mag,D2_difr,D2_mb]=calc_rates_magbins2(obs,Mmin,Mmax,dM);
[SR_forecast_magbin,observed_magbin,mag,D3_difr,D3_mb]=calc_rates_magbins2(SR_D3,Mmin,Mmax,dM);
[E5_forecast_magbin,observed_magbin,mag,D4_difr,D4_mb]=calc_rates_magbins2(E5_D3,Mmin,Mmax,dM);
[COM_D3_forecast_magbin,observed_magbin,mag,SR_difr,SR_mb]=calc_rates_magbins2(COM_D3,Mmin,Mmax,dM);


% Compare forecasted and observed seismicity for each magnitude bin
figure;
plot(mag(2:end),OBS_magbin,'r','LineWidth',8); hold on;
plot(mag(2:end),SR_forecast_magbin,'Color',[0.4 0.4 0.4],'LineWidth',7); hold on;grid on;
plot(mag(2:end),E5_forecast_magbin,'-','Color',[0.2 0.2 0.2],'LineWidth',2); hold on; 
plot(mag(2:end),COM_D3_forecast_magbin,'b','LineWidth',8); hold on;
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M_w','FontSize',14); 
% legend('obs','SR','E5','COMB');
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times')
set(get(gca,'YLabel'),'String','Rates within each magnitude bin','FontSize',18,'FontName','Times')
set(get(gca,'XLabel'),'String','M_w','FontSize',18,'FontName','Times')









