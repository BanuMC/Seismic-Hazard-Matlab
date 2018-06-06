clc
clear all
% close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% LOADING MODELS & OBSERVED DATA %%%%%%%%%%%%%%%%%%%%% 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Sum of Log likelihood for 15 days
ll_s0=-197.41;
ll_s1=-201.53;
ll_s2=-236.74;
ll_s3=-197.41;
ll_sr=-202.2;
ll_r0=-172.90;
ll_r1=-1431;
ll_r2=-635.86;
ll_e1=-426.4;
ll_e2=-346.7;
ll_e3=-220.14;
ll_e4=-233.83;
ll_e5=-249.77;


% Log likelihoods for 6 hours time bins
load loglike.mat;

log_like_S0=loglike.s0;
log_like_S1=loglike.s1;
log_like_S2=loglike.s2;
log_like_S3=loglike.s3;
log_like_SR=loglike.sr;
log_like_R0=loglike.r0;
log_like_R1=loglike.r1;
log_like_R2=loglike.r2;
log_like_E1=loglike.e1;
log_like_E2=loglike.e2;
log_like_E3=loglike.e3;
log_like_E4=loglike.e4;
log_like_E5=loglike.e5;




% figure;
% plot(log_like_SR);
% 
% 

log_like_E4(1)=-0.1;
log_like_E5(1)=-0.1;


like_SR=cumsum(log_like_SR);
like_R2=cumsum(log_like_R2);
like_E4=cumsum(log_like_E4);
like_E5=cumsum(log_like_E5);

%% CONSIDER SR, R2, E4 & E5 FOR AIC
% 15 days likelohoods
AIC_SR=2.*log(-1.*like_SR);
AIC_R2=2.*log(-1.*like_R2);
AIC_E4=2.*log(-1.*like_E4);
AIC_E5=2.*log(-1.*like_E5);


for i=1:length(log_like_SR);
    AIC_max(i)=max([AIC_SR(i) AIC_R2(i) AIC_E4(i) AIC_E5(i)]);
    sumAIC = 0;
    sumAIC = sumAIC + exp(-0.5*(AIC_SR(i)-AIC_max(i)));
    sumAIC = sumAIC + exp(-0.5*(AIC_R2(i)-AIC_max(i)));
    sumAIC = sumAIC + exp(-0.5*(AIC_E4(i)-AIC_max(i)));
    sumAIC = sumAIC + exp(-0.5*(AIC_E5(i)-AIC_max(i)));
    
    w1(i) = exp(-0.5*(AIC_SR(i)-AIC_max(i)))/sumAIC;
    w2(i) = exp(-0.5*(AIC_R2(i)-AIC_max(i)))/sumAIC;
    w3(i) = exp(-0.5*(AIC_E4(i)-AIC_max(i)))/sumAIC;
    w4(i) = exp(-0.5*(AIC_E5(i)-AIC_max(i)))/sumAIC;
    
    wg(i).a=[w1(i) w2(i) w3(i) w4(i)];
    clear sumAIC
end


A=zeros(length(log_like_SR),4);
A(1,:)=[0.25 0.25 0.25 0.25];
for j=2:length(log_like_SR);
    A(j,:)=wg(j-1).a;
end

AW1=A(:,1); % weight vector for SR model
AW2=A(:,2); % weight vector for R2 model
AW3=A(:,3); % weight vector for E4 model
AW4=A(:,4); % weight vector for E5 model

%% CALCULATE COMBINED MOEL
com = AW1.*SRR + AW2.*r22 + AW3.*e44 + AW4.*e55;

%% PLOTING
dt=0.25;
t=0.01:dt:length(log_like_SR)*dt;


figure;
% subplot(2,1,1);
plot(t,A(:,1),'r-o','LineWidth',2); hold on;
plot(t,A(:,2),'-o','Color',[0.4 0.4 0.4],'LineWidth',2); hold on;
plot(t,A(:,3),'k-o','LineWidth',2); hold on;
plot(t,A(:,4),'b-o','LineWidth',2); hold on;
grid on;
legend('SR','R2','E4','E5');
% xlabel('Time (days)','FontSize',14);
% ylabel('Akaike Weights','FontSize',14);
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',18,'FontName','Times')
set(get(gca,'YLabel'),'String','Akaike Weights','FontSize',18,'FontName','Times')


% figure;
% plot(t,SRR,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
% plot(t,r22,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
% plot(t,e44,'k','LineWidth',10); hold on;
% plot(t,e55,'b','LineWidth',3); hold on;
% plot(t,com,'g','LineWidth',4); hold on;
% plot(t,obs_rate,'r','LineWidth',4);  grid on;
% legend('SR','R2','E4','E5','Combined','observed');
% ylabel('Rates within each time bin','FontSize',14);
% xlabel('Time (days)','FontSize',14);



% save ('A.mat');

% save com.mat;







