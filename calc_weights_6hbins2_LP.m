clc
clear all
close all

load E5_LP.mat;
load SR_LP.mat;
load loglike_LP.mat;

% Observed rates
load rate_obs15d.mat;
obs_rate=rate_obs15d(1:60);


e5_lamda1=E5_LP(1).MODELS; % LP = 1 day
e5_lamda2=E5_LP(2).MODELS; % LP = 1.25 day
e5_lamda3=E5_LP(3).MODELS; % LP = 1.5 day
e5_lamda4=E5_LP(4).MODELS; % LP = 1.75 day
e5_lamda5=E5_LP(5).MODELS; % LP = 2 day
e5_lamda6=E5_LP(6).MODELS; % LP = 2.25 day
e5_lamda7=E5_LP(7).MODELS; % LP = 2.5 day
e5_lamda8=E5_LP(8).MODELS; % LP = 2.75 day
e5_lamda9=E5_LP(9).MODELS; % LP = 3 day
e5_lamda10=E5_LP(10).MODELS; % LP = 3.25 day
e5_lamda11=E5_LP(11).MODELS; % LP = 3.5 day
e5_lamda12=E5_LP(12).MODELS; % LP = 3.75 day
e5_lamda13=E5_LP(13).MODELS; % LP = 4.0 day


sr_lamda1=SR_LP(1).MODELS; % LP = 1 day
sr_lamda2=SR_LP(2).MODELS; % LP = 1.25 day
sr_lamda3=SR_LP(3).MODELS; % LP = 1.5 day
sr_lamda4=SR_LP(4).MODELS; % LP = 1.75 day
sr_lamda5=SR_LP(5).MODELS; % LP = 2 day
sr_lamda6=SR_LP(6).MODELS; % LP = 2.25 day
sr_lamda7=SR_LP(7).MODELS; % LP = 2.5 day
sr_lamda8=SR_LP(8).MODELS; % LP = 2.75 day
sr_lamda9=SR_LP(9).MODELS; % LP = 3 day
sr_lamda10=SR_LP(10).MODELS; % LP = 3.25 day
sr_lamda11=SR_LP(11).MODELS; % LP = 3.5 day
sr_lamda12=SR_LP(12).MODELS; % LP = 3.75 day
sr_lamda13=SR_LP(13).MODELS; % LP = 4.0 day



for ii=1:length(obs_rate);
    log_like_S1(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda1(ii))));
    log_like_S2(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda2(ii))));
    log_like_S3(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda3(ii))));
    log_like_S4(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda4(ii))));
    log_like_S5(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda5(ii))));
    log_like_S6(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda6(ii))));
    log_like_S7(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda7(ii))));
    log_like_S8(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda8(ii))));
    log_like_S9(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda9(ii))));
    log_like_S10(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda10(ii))));
    log_like_S11(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda11(ii))));
    log_like_S12(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda12(ii))));
    log_like_S13(ii) = sum(log(poisspdf(obs_rate(ii),sr_lamda13(ii))));
    log_like_E1(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda1(ii))));
    log_like_E2(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda2(ii))));
    log_like_E3(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda3(ii))));
    log_like_E4(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda4(ii))));
    log_like_E5(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda5(ii))));
    log_like_E6(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda6(ii))));
    log_like_E7(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda7(ii))));
    log_like_E8(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda8(ii))));
    log_like_E9(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda9(ii))));
    log_like_E10(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda10(ii))));
    log_like_E11(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda11(ii))));
    log_like_E12(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda12(ii))));
    log_like_E13(ii) = sum(log(poisspdf(obs_rate(ii),e5_lamda13(ii))));
end



[com1,t,log_like_com1,jll_com1,Ltest1]=calc_comb_LP(log_like_S1, log_like_E1, sr_lamda1, e5_lamda1,obs_rate);
[com2,t,log_like_com2,jll_com2,Ltest2]=calc_comb_LP(log_like_S2, log_like_E2, sr_lamda2, e5_lamda2,obs_rate);
[com3,t,log_like_com3,jll_com3,Ltest3]=calc_comb_LP(log_like_S3, log_like_E3, sr_lamda3, e5_lamda3,obs_rate);
[com4,t,log_like_com4,jll_com4,Ltest4]=calc_comb_LP(log_like_S4, log_like_E4, sr_lamda4, e5_lamda4,obs_rate);
[com5,t,log_like_com5,jll_com5,Ltest5]=calc_comb_LP(log_like_S5, log_like_E5, sr_lamda5, e5_lamda5,obs_rate);
[com6,t,log_like_com6,jll_com6,Ltest6]=calc_comb_LP(log_like_S6, log_like_E6, sr_lamda6, e5_lamda6,obs_rate);
[com7,t,log_like_com7,jll_com7,Ltest7]=calc_comb_LP(log_like_S7, log_like_E7, sr_lamda7, e5_lamda7,obs_rate);
[com8,t,log_like_com8,jll_com8,Ltest8]=calc_comb_LP(log_like_S8, log_like_E8, sr_lamda8, e5_lamda8,obs_rate);
[com9,t,log_like_com9,jll_com9,Ltest9]=calc_comb_LP(log_like_S9, log_like_E9, sr_lamda9, e5_lamda9,obs_rate);
[com10,t,log_like_com10,jll_com10,Ltest10]=calc_comb_LP(log_like_S10, log_like_E10, sr_lamda10, e5_lamda10,obs_rate);
[com11,t,log_like_com11,jll_com11,Ltest11]=calc_comb_LP(log_like_S11, log_like_E11, sr_lamda11, e5_lamda11,obs_rate);
[com12,t,log_like_com12,jll_com12,Ltest12]=calc_comb_LP(log_like_S12, log_like_E12, sr_lamda12, e5_lamda12,obs_rate);
[com13,t,log_like_com13,jll_com13,Ltest13]=calc_comb_LP(log_like_S13, log_like_E13, sr_lamda13, e5_lamda13,obs_rate);



y=[jll_com1 jll_com2 jll_com3 jll_com4 jll_com5 jll_com6 jll_com7 jll_com8 jll_com9 jll_com10 jll_com11 jll_com12 jll_com13];
x=1:13;
% figure;
% plot(x,y,'k'); hold on;
% plot(x,y,'ko');


ll_E5=loglike_LP.E5;
ll_SR=loglike_LP.SR;
ll_E5ref=loglike_LP.E5ref;
ll_SRref=loglike_LP.SRref;


t=1:0.25:4.0;
figure;

plot(t,ll_E5,'bo-','LineWidth',2); hold on;
plot(t,ll_SR,'ko-','LineWidth',2); hold on;
plot(t,y,'r^-','MarkerSize',8,'LineWidth',2);
plot(t,ll_E5ref,'k--','LineWidth',6); hold on;
plot(t,ll_SRref,'b--','LineWidth',6);
legend('E5','SR','COMB','E5-reference','SR-reference');
xlabel('Learning Period (days)','FontSize',18);
ylabel('log-likelihood','FontSize',18);
grid on;


ti=0.25:0.25:15;
figure;
plot(ti,e5_lamda9); hold on;
plot(ti,sr_lamda9,'k'); hold on;
plot(ti,com9,'g'); hold on;
plot(ti,obs_rate,'r');
legend('E5','SR','COMB','OBS');

% 
% COMB_MODEL_LP.com1=com1;
% COMB_MODEL_LP.com2=com2;
% COMB_MODEL_LP.com3=com3;
% COMB_MODEL_LP.com4=com4;
% COMB_MODEL_LP.com5=com5;
% COMB_MODEL_LP.com6=com6;
% COMB_MODEL_LP.com7=com7;
% COMB_MODEL_LP.com8=com8;
% COMB_MODEL_LP.com9=com9;
% COMB_MODEL_LP.com10=com10;
% COMB_MODEL_LP.com11=com11;
% COMB_MODEL_LP.com12=com12;
% COMB_MODEL_LP.com13=com13;
% savefile='COMB_MODEL_LP.mat';
% save(savefile,'COMB_MODEL_LP');




% loglike_LP_comb=y; %joint log likelihood of combined model for 1:0.25:4 days learning periods
% savefile='loglike_LP_comb.mat';
% save(savefile,'loglike_LP_comb');
% 
% Ltest_LP_comb=[Ltest1 Ltest2 Ltest3 Ltest4 Ltest5 Ltest6 Ltest7 Ltest8 Ltest9 Ltest10 Ltest11 Ltest12 Ltest13];
% savefile='Ltest_LP_comb.mat';
% save(savefile,'Ltest_LP_comb');


