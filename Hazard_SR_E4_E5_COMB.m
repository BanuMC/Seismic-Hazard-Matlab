clc
clear all
close all


%% Corinne's models
load all1907.mat
e4=rates_etas_fi_af_fb_09(:,1);
e5=allforecasts_e5_20d(:,3);
r2=rates_rt_tinc_09(:,1);
e44=e4(1:60);
e55=e5(1:60);
r22=r2(1:60);

% Observed rates
load rate_obs15d.mat;
obs_rate=rate_obs15d(1:60);


%% Shapiro Models: S0, S1, S2, S3, SR
load S0123R_15days_2.mat;  
SR=S0123R_15days_2.SR; SRR=SR(1:60); SRR=SRR';


savefile='SRR.mat';
save(savefile,'SRR');

load com.mat;
load A.mat;

time=0.25:0.25:(length(SRR))*0.25;


% figure;
% plot(time,SRR,'Color',[0.4 0.4 0.4],'LineWidth',7); hold on;
% plot(time,e44,'k','LineWidth',10); hold on;
% plot(time,e55,'b','LineWidth',3); hold on;
% plot(time,combo,'m-.','LineWidth',3); hold on;
% legend('SR','E4','E5','Combined');
% ylabel('Rates within 6 hours','FontSize',14);
% xlabel('Time (days)','FontSize',14);



% load b_upd.mat;
% b=b_upd(1:(length(SRR)));
% b=b';

b=1.5;

Mmin=0.9;
Mmax=5;
dm=0.1;




%% SR
mBack2=[];
bgmax2=[];
for m=Mmin:dm:Mmax;
    %if input magnitude changes, change here
    diff=0.9-m;
    bgm2=(SRR.*(10.^(b.*diff)));
    bgmax2=[bgmax2 bgm2];
end
k2=length(bgmax2(1,:))+2;
mBack2(:,3:k2)=bgmax2(:,:);
mBack2(:,1)=47.6;mBack2(:,2)= 7.5; %% arbitrarily selected location, corresponds to the r=0 distance
scpt_3=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack2(i,:),3,Mmax,Mmin,dm);
    scpt_3=[scpt_3 scP];
end
scpt_4=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack2(i,:),4,Mmax,Mmin,dm);
    scpt_4=[scpt_4 scP];
end
scpt_5=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack2(i,:),5,Mmax,Mmin,dm);
    scpt_5=[scpt_5 scP];
end


%% E4
mBack=[];
bgmax=[];
%%addjust line 25,26 in makeSwissbackgroundHazard_new for different min mag,max mag
for m=Mmin:dm:Mmax;
    %if input magnitude changes, change here
    diff=Mmin-m;
    bgm=(e44.*(10.^(b.*diff)));
    bgmax=[bgmax bgm];
end
k=length(bgmax(1,:))+2;
mBack(:,3:k)=bgmax(:,:);
mBack(:,1)=47.6;mBack(:,2)= 7.5;
scpint_3=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack(i,:),3,Mmax,Mmin,dm);
    scpint_3=[scpint_3 scP];
end
scpint_4=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack(i,:),4,Mmax,Mmin,dm);
    scpint_4=[scpint_4 scP];
end
scpint_5=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack(i,:),5,Mmax,Mmin,dm);
    scpint_5=[scpint_5 scP];
end


%% E5 
mBack3=[];
bgmax3=[];
for m=Mmin:dm:Mmax;
    %if input magnitude changes, change here
    diff=Mmin-m;
    bgm3=(e55.*(10.^(b.*diff)));
    bgmax3=[bgmax3 bgm3];
end
k3=length(bgmax3(1,:))+2;
mBack3(:,3:k2)=bgmax3(:,:);
mBack3(:,1)=47.6;mBack3(:,2)= 7.5; %% arbitrarily selected location, corresponds to the r=0 distance
sp_3=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),3,Mmax,Mmin,dm);
    sp_3=[sp_3 scP];
end
sp_4=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),4,Mmax,Mmin,dm);
    sp_4=[sp_4 scP];
end
sp_5=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),5,Mmax,Mmin,dm);
    sp_5=[sp_5 scP];
end





%% R2
mBack7=[];
bgmax7=[];
%%addjust line 25,26 in makeSwissbackgroundHazard_new for different min mag,max mag
for m=Mmin:dm:Mmax;
    %if input magnitude changes, change here
    diff=Mmin-m;
    bgm7=(r22.*(10.^(b.*diff)));
    bgmax7=[bgmax7 bgm7];
end
k7=length(bgmax7(1,:))+2;
mBack7(:,3:k)=bgmax7(:,:);
mBack7(:,1)=47.6;mBack7(:,2)= 7.5;
s_3=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP7,logP] = makeSwissBackgroundHazGrid_new(mBack7(i,:),3,Mmax,Mmin,dm);
    s_3=[s_3 scP7];
end
s_4=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP7,logP] = makeSwissBackgroundHazGrid_new(mBack7(i,:),4,Mmax,Mmin,dm);
    s_4=[s_4 scP7];
end
s_5=[];
for i=1:length(mBack(:,1))
    [scgx,scgy,scP7,logP] = makeSwissBackgroundHazGrid_new(mBack7(i,:),5,Mmax,Mmin,dm);
    s_5=[s_5 scP7];
end












%% COMBINED MODEL
mBack3=[];
bgmax3=[];
for m=Mmin:dm:Mmax;
    %if input magnitude changes, change here
    diff=Mmin-m;
    bgm3=(com.*(10.^(b.*diff)));
    bgmax3=[bgmax3 bgm3];
end
k3=length(bgmax3(1,:))+2;
mBack3(:,3:k2)=bgmax3(:,:);
mBack3(:,1)=47.6;mBack3(:,2)= 7.5; %% arbitrarily selected location, corresponds to the r=0 distance
sc_3=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),3,Mmax,Mmin,dm);
    sc_3=[sc_3 scP];
end
sc_4=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),4,Mmax,Mmin,dm);
    sc_4=[sc_4 scP];
end
sc_5=[];
for i=1:length(mBack2(:,1))
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_new(mBack3(i,:),5,Mmax,Mmin,dm);
    sc_5=[sc_5 scP];
end







% %% Ploting   
timel=0.25:0.25:length(SRR)*0.25;
% figure;
% % SR
% plot(timel,scpt_3,'k','LineWidth',2); hold on;
% plot(timel,scpt_4,'k--','LineWidth',2); hold on;
% plot(timel,scpt_5,'k.-','LineWidth',1); hold on;
% 
% % R2
% plot(timel,s_3,'r','LineWidth',2); hold on;
% plot(timel,s_4,'r--','LineWidth',2); hold on;
% plot(timel,s_5,'r.-','LineWidth',1); hold on;
% 
% % % E4
% % plot(timel,scpint_3,'r','LineWidth',2); hold on;
% % plot(timel,scpint_4,'r--','LineWidth',2); hold on;
% % plot(timel,scpint_5,'r.-','LineWidth',1); hold on;
% 
% % E5
% plot(timel,sp_3,'b','LineWidth',2); hold on;
% plot(timel,sp_4,'b--','LineWidth',2); hold on;
% plot(timel,sp_5,'b.-','LineWidth',1); hold on;
% 
% 
% plot(timel,sc_3,'Color',[0.4 0.4 0.4],'LineWidth',4); hold on;
% plot(timel,sc_4,'--','Color',[0.4 0.4 0.4],'LineWidth',4); hold on;
% plot(timel,sc_5,'.-','Color',[0.4 0.4 0.4],'LineWidth',4); hold on;
% 
% 
% 
% % legend('EMS3 - SR','EMS4 - SR','EMS5 - SR','EMS3 - R2','EMS4 - R2','EMS5 - R2','EMS3 - E5','EMS4 - E5','EMS5 - E5');
% legend('EMS3 - SR','EMS4 - SR','EMS5 - SR','EMS3 - R2','EMS4 - R2','EMS5 - R2','EMS3 - E5','EMS4 - E5','EMS5 - E5','EMS3 - COMB.','EMS4 - COMB.','EMS5 - COMB.');
% ylabel('Probability of exceeding EMS','FontSize',16);
% Xlabel('Time (days)','FontSize',16);
% % saveas(gcf, 'Prob_time_b1.5', 'jpg');



% %% STANDARD DEVIATION AND 95% CI
% sigma_3=std(sc_3); % standard deviation
% CI_3=(2*sigma_3)/(sqrt(length(sc_3))); % 95% CI
% 
% sigma_4=std(sc_4); % standard deviation
% CI_4=(2*sigma_4)/(sqrt(length(sc_4))); % 95% CI
% 
% sigma_5=std(sc_5); % standard deviation
% CI_5=(2*sigma_5)/(sqrt(length(sc_5))); % 95% CI
% 
% a=sc_3+CI_3;
% b=sc_3-CI_3


% figure;
% plot(timel,sc_3,'k','LineWidth',2); hold on;
% % plot(timel,sc_3+CI_3,'b','LineWidth',1); hold on;
% % plot(timel,sc_3-CI_3,'b','LineWidth',1); hold on;
% plot(timel,sc_4,'k--','LineWidth',2); hold on;
% plot(timel,sc_5,'k.-','LineWidth',1); hold on;
% 
% % ylabel('Probability of exceeding EMS','FontSize',16);
% % Xlabel('Time (days)','FontSize',16);
% legend('EMS3','EMS4','EMS5');





w1=A(:,1);
w2=A(:,2);
w3=A(:,3);
w4=A(:,4);


comb_haz_ems3 = (scpt_3'.*w1) + (s_3'.*w2) + (scpint_3'.*w3) + (sp_3'.*w4);
comb_haz_ems4 = (scpt_4'.*w1) + (s_4'.*w2) + (scpint_4'.*w3) + (sp_4'.*w4);
comb_haz_ems5 = (scpt_5'.*w1) + (s_5'.*w2) + (scpint_5'.*w3) + (sp_5'.*w4);

% figure;
% plot(timel,sc_3,'k','LineWidth',2); hold on;
% plot(timel,comb_haz,'r','LineWidth',2); 


figure;
subplot(2,1,1);
plot(timel,scpt_4,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
plot(timel,s_4,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(timel,scpint_4,'Color',[0.8 0.8 0.8],'LineWidth',8); hold on;
plot(timel,sp_4,'--','Color',[0.1 0.1 0.1],'LineWidth',2); hold on;
plot(timel,comb_haz_ems4,'k','LineWidth',6); 
legend('SR','R2','E4','E5','COMB');
ylim([0 1]);
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'YLabel'),'String','Probability of exceeding EMS V','FontSize',24,'FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',24,'FontName','Times')


subplot(2,1,2);
plot(timel,scpt_3,'Color',[0.4 0.4 0.4],'LineWidth',14); hold on;
plot(timel,s_3,'Color',[0.2 0.2 0.2],'LineWidth',2); hold on;
plot(timel,scpint_3,'Color',[0.8 0.8 0.8],'LineWidth',8); hold on;
plot(timel,sp_3,'Color',[0.1 0.1 0.1],'LineWidth',2); hold on;
plot(timel,comb_haz_ems3,'k','LineWidth',6); 
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'YLabel'),'String','Probability of exceeding EMS IV','FontSize',24,'FontName','Times')
set(get(gca,'XLabel'),'String','Time [days]','FontSize',24,'FontName','Times')






% figure;
% % subplot(1,3,3);
% plot(timel,scpt_5,'k','LineWidth',2); hold on;
% plot(timel,s_5,'r','LineWidth',2); hold on;
% plot(timel,scpint_5,'m','LineWidth',2); hold on;
% plot(timel,sp_5,'b','LineWidth',2); hold on;
% plot(timel,comb_haz_ems5,'Color',[0.4 0.4 0.4],'LineWidth',2); 
% ylabel('Probability of exceeding EMS V','FontSize',16);
% Xlabel('Time (days)','FontSize',16);
% grid on;
% ylim([0 1]);
% 
% legend('SR','R2','E4','E5','COMB');




