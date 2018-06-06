clc
clear all
close all

%% Corinne's models
load all1907.mat
r0=rates_copost(:,1);
r1=rates_gen2(:,1);
r2=rates_rt_tinc_09(:,1);

e1=allforecasts_fixacpb_025(:,3);
e2=allforecasts_fixapcb_fi_20d(:,3);
e3=allforecasts_fixb_cst_20d(:,3);
e4=rates_etas_fi_af_fb_09(:,1);
e5=allforecasts_e5_20d(:,3);

r00=r0(1:60);
r11=r1(1:60);
r22=r2(1:60);
e11=e1(1:60);
e22=e2(1:60);
e33=e3(1:60);
e44=e4(1:60);
e55=e5(1:60);

% Observed rates
load rate_obs15d.mat;
obs_rate=rate_obs15d(1:60);


%% Shapiro Models: S0, S1, S2, S3, SR
load S0123R_15days_2.mat;  

S0=S0123R_15days_2.S0; SS0=S0(1:60); SS0=SS0';
S1=S0123R_15days_2.S1; SS1=S1(1:60); SS1=SS1';
S2=S0123R_15days_2.S2; SS2=S2(1:60); SS2=SS2'; 
S3=S0123R_15days_2.S3; SS3=S3(1:60); SS3=SS3';
SR=S0123R_15days_2.SR; SRR=SR(1:60); SRR=SRR';

load com.mat;



Mmin=0.9;
Mmax=3.5;
dM=0.1;

%% Calculate rates for magnitude bins
% S MODELS
[S0_forecast_magbin,observed_magbin,mag,S0_difr,S0_mb]=calc_rates_magbins2(SS0,Mmin,Mmax,dM);
[S1_forecast_magbin,observed_magbin,mag,S1_difr,S1_mb]=calc_rates_magbins2(SS1,Mmin,Mmax,dM);
[S2_forecast_magbin,observed_magbin,mag,S2_difr,S2_mb]=calc_rates_magbins2(SS2,Mmin,Mmax,dM);
[S3_forecast_magbin,observed_magbin,mag,S3_difr,S3_mb]=calc_rates_magbins2(SS3,Mmin,Mmax,dM);
[SR_forecast_magbin,observed_magbin,mag,SR_difr,SR_mb]=calc_rates_magbins2(SRR,Mmin,Mmax,dM);
% R MODELS
[R0_forecast_magbin,observed_magbin,mag,R0_difr,R0_mb]=calc_rates_magbins2(r00,Mmin,Mmax,dM);
[R1_forecast_magbin,observed_magbin,mag,R1_difr,R1_mb]=calc_rates_magbins2(r11,Mmin,Mmax,dM);
[R2_forecast_magbin,observed_magbin,mag,R2_difr,R2_mb]=calc_rates_magbins2(r22,Mmin,Mmax,dM);
% E MODELS
[E1_forecast_magbin,observed_magbin,mag,E1_difr,E1_mb]=calc_rates_magbins2(e11,Mmin,Mmax,dM);
[E2_forecast_magbin,observed_magbin,mag,E2_difr,E2_mb]=calc_rates_magbins2(e22,Mmin,Mmax,dM);
[E3_forecast_magbin,observed_magbin,mag,E3_difr,E3_mb]=calc_rates_magbins2(e33,Mmin,Mmax,dM);
[E4_forecast_magbin,observed_magbin,mag,E4_difr,E4_mb]=calc_rates_magbins2(e44,Mmin,Mmax,dM);
[E5_forecast_magbin,observed_magbin,mag,E5_difr,E5_mb]=calc_rates_magbins2(e55,Mmin,Mmax,dM);
% COMBO
[C_forecast_magbin,observed_magbin,mag,COMBO_difr,COMBO_mb]=calc_rates_magbins2(com,Mmin,Mmax,dM);

%% NUMBER OF EVENTS WITH M>=2
mb2_obs=sum(observed_magbin(12:21))
% S MODELS
mb2_S0=sum(S0_forecast_magbin(12:21))
mb2_S1=sum(S1_forecast_magbin(12:21))
mb2_S2=sum(S2_forecast_magbin(12:21))
mb2_S3=sum(S3_forecast_magbin(12:21))
mb2_SR=sum(SR_forecast_magbin(12:21))
% R MODELS
mb2_R0=sum(R0_forecast_magbin(12:21))
mb2_R1=sum(R1_forecast_magbin(12:21))
mb2_R2=sum(R2_forecast_magbin(12:21))
% E MODELs
mb2_E1=sum(E1_forecast_magbin(12:21))
mb2_E2=sum(E2_forecast_magbin(12:21))
mb2_E3=sum(E3_forecast_magbin(12:21))
mb2_E4=sum(E4_forecast_magbin(12:21))
mb2_E5=sum(E5_forecast_magbin(12:21))



% %% NUMBER OF EVENTS WITH [0.9 <= M < 1.5]
% mbin1_obs=sum(observed_magbin(1:6))
% % S MODELS
% mbin1_S0=sum(S0_forecast_magbin(1:6))-mbin1_obs
% mbin1_S1=sum(S1_forecast_magbin(1:6))-mbin1_obs
% mbin1_S2=sum(S2_forecast_magbin(1:6))-mbin1_obs
% mbin1_S3=sum(S3_forecast_magbin(1:6))-mbin1_obs
% mbin1_SR=sum(SR_forecast_magbin(1:6))-mbin1_obs
% % R MODELS
% mbin1_R0=sum(R0_forecast_magbin(1:6))-mbin1_obs
% mbin1_R1=sum(R1_forecast_magbin(1:6))-mbin1_obs
% mbin1_R2=sum(R2_forecast_magbin(1:6))-mbin1_obs
% % E MODELS
% mbin1_E1=sum(E1_forecast_magbin(1:6))-mbin1_obs
% mbin1_E2=sum(E2_forecast_magbin(1:6))-mbin1_obs
% mbin1_E3=sum(E3_forecast_magbin(1:6))-mbin1_obs
% mbin1_E4=sum(E4_forecast_magbin(1:6))-mbin1_obs
% mbin1_E5=sum(E5_forecast_magbin(1:6))-mbin1_obs



% %% NUMBER OF EVENTS WITH [1.5 <= M < 2]
% mbin2_obs=sum(observed_magbin(7:11))
% % S MODELS
% mbin2_S0=sum(S0_forecast_magbin(7:11))-mbin2_obs
% mbin2_S1=sum(S1_forecast_magbin(7:11))-mbin2_obs
% mbin2_S2=sum(S2_forecast_magbin(7:11))-mbin2_obs
% mbin2_S3=sum(S3_forecast_magbin(7:11))-mbin2_obs
% mbin2_SR=sum(SR_forecast_magbin(7:11))-mbin2_obs
% % R MODELS
% mbin2_R0=sum(R0_forecast_magbin(7:11))-mbin2_obs
% mbin2_R1=sum(R1_forecast_magbin(7:11))-mbin2_obs
% mbin2_R2=sum(R2_forecast_magbin(7:11))-mbin2_obs
% % E MODELS
% mbin2_E1=sum(E1_forecast_magbin(7:11))-mbin2_obs
% mbin2_E2=sum(E2_forecast_magbin(7:11))-mbin2_obs
% mbin2_E3=sum(E3_forecast_magbin(7:11))-mbin2_obs
% mbin2_E4=sum(E4_forecast_magbin(7:11))-mbin2_obs
% mbin2_E5=sum(E5_forecast_magbin(7:11))-mbin2_obs


%% COMPUTE LOG-LIKELIHOODS
for ii=1:length(observed_magbin);
    % S MODELS
    log_S0(ii) = sum(log(poisspdf(observed_magbin(ii),S0_forecast_magbin(ii))));
    log_S1(ii) = sum(log(poisspdf(observed_magbin(ii),S1_forecast_magbin(ii))));
    log_S2(ii) = sum(log(poisspdf(observed_magbin(ii),S2_forecast_magbin(ii))));
    log_S3(ii) = sum(log(poisspdf(observed_magbin(ii),S3_forecast_magbin(ii))));
    log_SR(ii) = sum(log(poisspdf(observed_magbin(ii),SR_forecast_magbin(ii))));
    % R MODELS
    log_R0(ii) = sum(log(poisspdf(observed_magbin(ii),R0_forecast_magbin(ii))));
    log_R1(ii) = sum(log(poisspdf(observed_magbin(ii),R1_forecast_magbin(ii))));
    log_R2(ii) = sum(log(poisspdf(observed_magbin(ii),R2_forecast_magbin(ii))));
    % E MODELS
    log_E1(ii) = sum(log(poisspdf(observed_magbin(ii),E1_forecast_magbin(ii))));
    log_E2(ii) = sum(log(poisspdf(observed_magbin(ii),E2_forecast_magbin(ii))));
    log_E3(ii) = sum(log(poisspdf(observed_magbin(ii),E3_forecast_magbin(ii))));
    log_E4(ii) = sum(log(poisspdf(observed_magbin(ii),E4_forecast_magbin(ii))));
    log_E5(ii) = sum(log(poisspdf(observed_magbin(ii),E5_forecast_magbin(ii))));
    % COMBINED MODEL
    log_C(ii) = sum(log(poisspdf(observed_magbin(ii),C_forecast_magbin(ii))));
end

%% SUM OF LOG-LIKELIHOODS FOR 15 DAYS FORECASTS WITH SHAPIRO & ETAS & STEP MODELS
% S MODELS 
lls_S0=sum(log_S0)
lls_S1=sum(log_S1)
lls_S2=sum(log_S2)
lls_S3=sum(log_S3)
lls_SR=sum(log_SR)
% R MODELS
lls_R0=sum(log_R0)
lls_R1=sum(log_R1)
lls_R2=sum(log_R2)
% E MODELS
lls_E1=sum(log_E1)
lls_E2=sum(log_E2)
lls_E3=sum(log_E3)
lls_E4=sum(log_E4)
lls_E5=sum(log_E5)
% COMBINED MODEL
lls_C=sum(log_C)

%% PLOTING

time=0.25:0.25:(length(obs_rate))*0.25;

%% PLOT SHAPIRO MODELS
% %FORECASTS: S0, S1, S2, S3, SR
% figure; 
% plot(time,obs_rate,'r','LineWidth',4); hold on; grid on;
% plot(time,SS0,'Color',[0.4 0.4 0.4],'LineWidth',7);  hold on;
% plot(time,SS1,'-.','Color',[0.7 0.7 0.7],'LineWidth',9); hold on;
% plot(time,SS2,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(time,SS3,'k','LineWidth',2); hold on;
% plot(time,SRR,'b','LineWidth',2); 
% legend('obs','S0','S1','S2','S3','SR');
% ylabel('Rates within 6 hours','FontSize',14);
% xlabel('Time (days)','FontSize',14);
% % title('Shapiro Models');

% % Cumulative distribution following Gutenberg Richter
% figure;
% plot(mag,S0_mb); hold on; grid on; 
% plot(mag,S1_mb,'r'); hold on; 
% plot(mag,S2_mb,'c'); hold on;  
% plot(mag,S3_mb,'m'); hold on; 
% plot(mag,SR_mb,'g'); hold on; 
% xlabel('M');
% ylabel('Number of events N (N>M)');
% title('Magnitude Frequency distribution (Gutenberg Richter law)');
% legend('S0','S1','S2','S3','SR');

% % Compare forecasted and observed seismicity for each magnitude bin
% figure;
% plot(mag(2:end),observed_magbin,'r','LineWidth',4); hold on; grid on;
% plot(mag(2:end),S0_forecast_magbin,'Color',[0.4 0.4 0.4],'LineWidth',7); hold on; 
% plot(mag(2:end),S1_forecast_magbin,'-.','Color',[0.7 0.7 0.7],'LineWidth',9); hold on;
% plot(mag(2:end),S2_forecast_magbin,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(mag(2:end),S3_forecast_magbin,'k','LineWidth',2); hold on;
% plot(mag(2:end),SR_forecast_magbin,'b','LineWidth',2); hold on;
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M','FontSize',14); 
% legend('observed','S0','S1','S2','S3','SR');
% % title('Comparison of Shapiro Models','FontSize',14);

% % Plot the difference between observed seismicity and forecasted seismicity for each magnitude bin
% figure;
% plot(mag(2:end),S0_difr,'b-o'); hold on; grid on;
% plot(mag(2:end),S1_difr,'r-o'); hold on;
% plot(mag(2:end),S2_difr,'c-o'); hold on;
% plot(mag(2:end),S3_difr,'m-o'); hold on;
% plot(mag(2:end),SR_difr,'g-o'); hold on;
% xlabel('M'); 
% ylabel('N_f_o_r_e_c_a_s_t_e_d-N_o_b_s_e_r_v_e_d');
% legend('S0','S1','S2','S3','SR');

% %% PLOT R-MODELS' 
% % FORECASTS: R0, R1, R2 
% figure; 
% plot(time,obs_rate,'r','LineWidth',4); hold on; grid on;
% plot(time,r00,'Color',[0.7 0.7 0.7],'LineWidth',7);  hold on;
% plot(time,r11,'-.','Color',[0.4 0.4 0.4],'LineWidth',7); hold on;
% plot(time,r22,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% legend('obs','R0','R1','R2');
% ylabel('Rates within 6 hours','FontSize',14);
% xlabel('Time (days)','FontSize',14);
% 
% % Compare forecasted and observed seismicity for each magnitude bin
% figure;
% plot(mag(2:end),observed_magbin,'r','LineWidth',4); hold on; grid on;
% plot(mag(2:end),R0_forecast_magbin,'Color',[0.7 0.7 0.7],'LineWidth',7); hold on; 
% plot(mag(2:end),R1_forecast_magbin,'-.','Color',[0.4 0.4 0.4],'LineWidth',7); hold on;
% plot(mag(2:end),R2_forecast_magbin,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M','FontSize',14); 
% legend('observed','R0','R1','R2');


% %% PLOT E-MODELS
% %FORECASTS: E1, E2, E3, E4, E5
% figure; 
% plot(time,obs_rate,'r','LineWidth',4); hold on; grid on;
% plot(time,e11,'Color',[0.4 0.4 0.4],'LineWidth',7);  hold on;
% plot(time,e22,'-.','Color',[0.7 0.7 0.7],'LineWidth',9); hold on;
% plot(time,e33,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(time,e44,'k','LineWidth',2); hold on;
% plot(time,e55,'b','LineWidth',2); 
% legend('obs','E1','E2','E3','E4','E5');
% ylabel('Rates within 6 hours','FontSize',14);
% xlabel('Time (days)','FontSize',14);
% 
% % Compare forecasted and observed seismicity for each magnitude bin
% figure;
% plot(mag(2:end),observed_magbin,'r','LineWidth',4); hold on; grid on;
% plot(mag(2:end),E1_forecast_magbin,'Color',[0.4 0.4 0.4],'LineWidth',7); hold on; 
% plot(mag(2:end),E2_forecast_magbin,'-.','Color',[0.7 0.7 0.7],'LineWidth',9); hold on;
% plot(mag(2:end),E3_forecast_magbin,'-.','Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(mag(2:end),E4_forecast_magbin,'k','LineWidth',2); hold on;
% plot(mag(2:end),E5_forecast_magbin,'b','LineWidth',2); hold on;
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M','FontSize',14); 
% legend('observed','E1','E2','E3','E4','E5');





%% SELECTED MODELS: SR, R2, E4, E5
% Selection criterion:
% 1) Parameters updated 
% 2) Flowrate or injected fluid volume is included in the computations
figure;
subplot(2,1,1);
plot(time,SRR,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
plot(time,r22,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(time,e44,'b','LineWidth',6); hold on;
plot(time,e55,'c','LineWidth',3); hold on;
plot(time,com,'r','LineWidth',6); hold on;
plot(time,obs_rate,'k','LineWidth',6);  
legend('SR','R2','E4','E5','COMB','Obs');
ylim([0 300]);
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',20,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within 6 hours','FontSize',20,'FontName','Times')

subplot(2,1,2);
plot(mag(2:end),SR_forecast_magbin,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
plot(mag(2:end),R2_forecast_magbin,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(mag(2:end),E4_forecast_magbin,'b','LineWidth',6); hold on;
plot(mag(2:end),E5_forecast_magbin,'c','LineWidth',3); hold on;
plot(mag(2:end),C_forecast_magbin,'r','LineWidth',6); hold on;
plot(mag(2:end),observed_magbin,'k','LineWidth',4); hold on; 
% legend('SR','R2','E4','E5','Combined','observed');
% ylabel('Rates within each magnitude bin','FontSize',14);
% xlabel('M','FontSize',14); 
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','M_w','FontSize',20,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within each magnitude bin','FontSize',20,'FontName','Times')



all_models.SR=SRR;
all_models.R2=r22;
all_models.E4=e44;
all_models.E5=e55;
all_models.COMB=com;
all_models.obs=obs_rate;

savefile='all_models.mat';
save(savefile,'all_models');








% %% SELECTED MODELS: SR, R2, E4, E5
% % Selection criterion:
% % 1) Parameters updated 
% % 2) Flowrate or injected fluid volume is included in the computations
% figure;
% % subplot(2,1,1);
% plot(time,SRR,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
% plot(time,r22,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(time,e44,'k','LineWidth',10); hold on;
% plot(time,e55,'LineWidth',5,'Color',[.8 .8 .8]); hold on;
% plot(time,com,'LineWidth',9,'Color',[.5 .5 .5]); hold on;
% plot(time,obs_rate,'r','LineWidth',4);  grid on;
% legend('SR','R2','E4','E5','COMB','Obs');
% % ylabel('Rates within 6 hours','FontSize',14);
% % xlabel('Time (days)','FontSize',14);
% set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','Time [days]','FontSize',18,'FontName','Times')
% set(get(gca,'YLabel'),'String','Rates within 6 hours','FontSize',18,'FontName','Times')

% subplot(2,1,2);
% plot(mag(2:end),SR_forecast_magbin,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
% plot(mag(2:end),R2_forecast_magbin,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(mag(2:end),E4_forecast_magbin,'k','LineWidth',10); hold on;
% plot(mag(2:end),E5_forecast_magbin,'g','LineWidth',3); hold on;
% plot(mag(2:end),C_forecast_magbin,'y','LineWidth',6); hold on;
% plot(mag(2:end),observed_magbin,'r','LineWidth',4); hold on; grid on;
% % legend('SR','R2','E4','E5','Combined','observed');
% % ylabel('Rates within each magnitude bin','FontSize',14);
% % xlabel('M','FontSize',14); 
% set(gca,'LineWidth',2,'FontSize',18,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','M_w','FontSize',18,'FontName','Times')
% set(get(gca,'YLabel'),'String','Rates within each magnitude bin','FontSize',18,'FontName','Times')
% 







% figure;
% semilogy(mag(2:end),observed_magbin,'k','LineWidth',4); hold on; grid on;














% load 'basel_spimed_25708.mat';
% dstart=0.01;
% dstep=0.25;
% a(:,10)=0;
% 
% forecast_file_SRR=ones(60,10); 
% forecast_file_SRR(:,1)=7.5904;
% forecast_file_SRR(:,2)=7.6014;
% forecast_file_SRR(:,3)=47.5786;
% forecast_file_SRR(:,4)=47.5896;
% forecast_file_SRR(:,5)=0;
% forecast_file_SRR(:,6)=1;
% forecast_file_SRR(:,7)=0.9;
% forecast_file_SRR(:,8)=3;
% forecast_file_SRR(:,9)=SRR;
% forecast_file_SRR(:,10)=1;
% 
% [rResult_SRR]=NLtest_wrapper_new(a, forecast_file_SRR,dstart,dstep);
% Ltest_SRR=rResult_SRR.fNDailyRatioReject
% up=rResult_SRR.vLL+rResult_SRR.vLL_q975;
% down=rResult_SRR.vLL-rResult_SRR.vLL_q025;
% ll=rResult_SRR.fLL
% vll=rResult_SRR.vLL
% 
% figure;
% errorbar(time,vll,rResult_SRR.vLL_q975,rResult_SRR.vLL_q025); hold on; plot(time,ll,'ro');
% plot(time,ll,'ro'); hold on;
% % plot(time,up,'bo'); hold on;
% % plot(time,down,'go');
% xlabel('Time');
% ylabel('Log likelihood');
% % title('SR Model')
% 
% 
% forecast_file_E5(:,:)=forecast_file_SRR(:,:);
% forecast_file_E5(:,9)=e55;
% [rResult_E5]=NLtest_wrapper_new(a, forecast_file_E5,  dstart,dstep);
% Ltest_E5=rResult_E5.fNDailyRatioReject
% up_e5=rResult_E5.vLL+rResult_E5.vLL_q975;
% down_e5=rResult_E5.vLL-rResult_E5.vLL_q025;
% ll_e5=rResult_E5.fLL
% vll_e5=rResult_E5.vLL
% 
% figure;
% errorbar(time,vll_e5,rResult_E5.vLL_q975,rResult_E5.vLL_q025); hold on; plot(time,ll_e5,'ro');
% plot(time,ll_e5,'ro'); hold on;
% % plot(time,up_e5,'bo'); hold on;
% % plot(time,down_e5,'go');
% xlabel('Time');
% ylabel('Log likelihood');
% % title('E5 Model');
% 
% % models.RO=r00;
% % models.R1=r11;
% % models.R2=r22;
% % models.E1=e11;
% % models.E2=e22;
% % models.E3=e33;
% % models.E4=e44;
% % models.E5=e55;
% % models.obs=obs_rate;
% % 
% % savefile='models.mat';
% % save(savefile,'models');
% 





