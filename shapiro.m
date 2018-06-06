clc
clear all
close all

%% Read the catalog
fid = fopen(['/Users/banumenacabrera/GEISER/ETAS/DATA/cat.dat'],'r');
line=fgetl(fid);
data = fscanf(fid,'%e %e %e %e %e',[5 inf]);  
fclose(fid);
data = data';
%no = data(:,1); 
ctime = data(:,2); % catalog time
%lat = data(:,3); 
%lon = data(:,4);
mag = data(:,5); % catalog magnitude
%% Read the flow rate
fid2 = fopen(['/Users/banumenacabrera/GEISER/ETAS/DATA/fluxrate.dat'],'r');
line=fgetl(fid2);
data2 = fscanf(fid2,'%e %e %e %e %e %e %e %e',[8 inf]);  
fclose(fid2);
data2 = data2';
time = data2(:,2); % time in days
flow_rate = data2(:,7); % flow rate in liter/sec



%% INPUT
fForecastStart = 0.01; % Start of the forecast window
fForecastEnd = 15; % End of the forecast window
fDt     = 0.25; % Time step in days corresponding to 6 hours
% fDt     = 0.125; % Time step in days corresponding to 3 hours
% fDt     = 0.125/2; % Time step in days corresponding to 1.5 hours
Mmin=0.9;

Mcat=0.9;
% Compute b value from Basel data 
[b, b_upd, tb] =  calc_bvalue(Mcat,fForecastEnd,fDt); 

% Magn=mag';
% [fMeanMag1, fBValue1, fStdDev1, fAValue1] =  calc_bmemagMag(Magn);

save b_upd.mat;
% mean_b_upd=mean(b_upd)

% % Calculate injected fluid volume from flow rate
% % Calculated only once because computation takes too much time 
% tt(1)=0;
% for ii=1:(length(time)-1);
%     tt(ii+1)=time(ii+1)-time(ii);
% end
% flow_rate=flow_rate./10^3; % convert from liter to m3 (m3/min)
% tt=tt.*(24*60); % convert time [days] into time [minutes]
% fluid_volume = tt'.*flow_rate; % tt[minutes] * flowrate[m3/min] = injected fluid volume [m3]
% fluid_volume_total=sum(fluid_volume); % total injected fluid volume in m3
% fluid_volume_cumulative = cumsum(fluid_volume); % cumulative sum of fluid volume
% fvc=fluid_volume_cumulative;
% save fvc

%% Load injected fluid volume
fluid = load('fvc.mat');
fluid_volume_cumulative=fluid.fluid_volume_cumulative;

% The first 15 days of the catalog
t15=ctime(ctime<=15);
mag15=mag(1:length(t15));



m0=10.^((mag15+10.7).*1.5); % seismic moment in dyne cm (Hanks and Kanamori, 1979)
for ii=1:15;
    x=ii-t15;
    ind=find(x==min(x(x>0)));
    if length(ind)>1;
        ind=ind(end);
    end
    indic(ii)=ind;
    m(ii)=sum(m0(1:indic(ii)));
    clear x, ind;
end
timet=[1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15];
mm=[m(1) m(1) m(2) m(2) m(3) m(3) m(4) m(4) m(5) m(5) m(6) m(6) m(7) m(7) m(8) m(8) m(9) m(9) m(10) m(10) m(11) m(11) m(12) m(12) m(13) m(13) m(14) m(14) m(15)];








% Events with magnitudes greated than or equal to the catalog magnitude
MM=mag15(mag15>=Mcat); % magnitudes larger than catalog magnitude
ind=find(mag15>=Mcat); % indices of magnitudes >= Mcat
% Time vector corresponding to events of M>=Mcat
for i=1:length(ind);
    t(i)=t15(ind(i));
end
t=t';

magnitude=MM';
[fMeanMag, fBValue, fStdDev, fAValue] =  calc_bmemagMag(magnitude);

%% Compute observed seismicity rates from the catalog within 6 hours time intervals
k=1;
for ij=fForecastStart:fDt:fForecastEnd+fDt;
    fr(k)=ij+fDt;
    k=k+1;
end
rate_obs=[];
for ii=1:length(fr);
    sr_d = t <= fr(ii) & t > fr(ii)-fDt;
    rob=sum(sr_d==1);
    rate_obs=[rate_obs rob];
end
rate_obs=rate_obs'; % Observed seismicity rate within every 6 hours time bins 

% rate_obs15d=rate_obs;
% save rate_obs15d.mat;

%% Injected fluid volume within 6 hours time intervals
for ih=1:(12/fDt); %% we have no data after 12 days for injected fluid volume
    indic = find(time <= fr(ih) & time > fr(ih)-fDt);
    ifv(ih)=fluid_volume_cumulative(max(indic));
    clear indic;
end
len=length(ifv);
for iy=1:length(fr)-length(ifv);
    ifv(len+iy)=ifv(len);
end
ifvn=ifv(1:end-1);
delta_ifv(1)=ifvn(1);
for jk=1:length(ifvn)-1;
    delta_ifv(jk+1)=ifvn(jk+1)-ifvn(jk); % injected fluid volume within every 6 hours bins during 15 days
end
index=find(delta_ifv==0);
index1=min(index);

savefile='delta_ifv.mat';
save(savefile,'delta_ifv');



%% FOR UPDATING THE FORECAST, UPDATE THE SEISMOGENIC INDEX FOR 6, 12, 18, 24, ....600 HOURS
rate_obs_update=[];
% Number of events for every 6, 12, 18 24, ... 600 hours (0.25, 0.5, 0.75,..... 15 days) (rate_obs_update)
for ii=1:(length(fr)-1);
    sr_n = t <= fr(ii+1);
    rob_n=sum(sr_n==1);
    rate_obs_update=[rate_obs_update rob_n];
end

rate_obs_upd=rate_obs_update;
savefile='rate_obs_upd.mat';
save(savefile,'rate_obs_upd');




% Injected fluid volume for every 6, 12, 18 24, ... 600 hours (0.25, 0.5, 0.75,..... 15 days) (ifv2)
for ih=1:(12/fDt);
    indic2 = find(time <= fr(ih));
    ifv2(ih)=fluid_volume_cumulative(max(indic2));
    clear indic2
end
len2=length(ifv2);
for ig=1:length(fr)-length(ifv2)-1;
    ifv2(len2+ig)=ifv2(len2);
end

in_fv2=ifv2;
savefile='in_fv2.mat';
save(savefile,'in_fv2');

% Update seismogenic index every 6, 12, 18, 24, ....hours
eps_upd(1)=0.5;
for ip=2:length(ifv2);
    eps_upd(ip)=log10(rate_obs_update(ip))-log10(ifv2(ip))+b*Mcat; %SEISMOGENIC INDEX UPDATED AT 6, 12, 18 ... HOURS
end

figure; plot(ifv2);

%% PROBLEM, WHEN INJECTED FLUID VOLUME BECOMES ZERO, HOW TO MODEL LAMDA????

%% S1 MODEL
% SEISMOGENIC INDEX : FREE
% B VALUE : FIXED
% Forecast using updated seismogenic index, updated at 6, 12, 18, ... hours
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
         delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); %%%%%%%%%%%%%%%%% bu dogru degil, sadece deneme icin
    end
    lamda_s1(ic)=10^(log10(delta_ifv(ic))-(b*Mmin)+eps_upd(ic));
end

tbv=zeros(length(b_upd)+1,1);
bv=zeros(length(b_upd)+1,1);
tbv(1)=0.001;
tbv(2:end)=tb;
bv(1)=b_upd(1);
bv(2:end)=b_upd;

% Update seismogenic index every 6, 12, 18, 24, ....hours 
% Using updated b value
eps_upd2(1)=0.5;
for ip=2:length(ifv2);
    eps_upd2(ip)=log10(rate_obs_update(ip))-log10(ifv2(ip))+bv(ip)*Mcat; %SEISMOGENIC INDEX UPDATED AT 6, 12, 18 ... HOURS
end











%% SR MODEL ------> UPDATING b value as more data arrives 
% SEISMOGENIC INDEX : FREE
% B VALUE : FREE
% Forecast using updated seismogenic index, updated at 6, 12, 18, ... hours
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
        delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); %%%%%%%%%%%%%%%%% 
    end
    lamda_sr(ic)=10^(log10(delta_ifv(ic))-(bv(ic)*Mmin)+eps_upd2(ic));
end



savefile='lamda_sr.mat';
save(savefile,'lamda_sr');

in_fv=delta_ifv;

savefile='in_fv.mat';
save(savefile,'in_fv');


savefile='bv.mat';
save(savefile,'bv');

epsilon=eps_upd2;
savefile='epsilon.mat';
save(savefile,'epsilon');


M=2;
for ic=1:length(delta_ifv);
    P(ic)=exp(-1*delta_ifv(ic)*10^eps_upd2(ic)-bv(ic)*M);
end



%% S0 MODEL
% SEISMOGENIC INDEX : FIXED
% B VALUE : FIXED (BASED ON BASEL DATA)
eps_global=(log10(rate_obs_update(23)))-(log10(ifv2(23)))+(b*Mcat)
% Forecast using global seismogenic index
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
         delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); 
    end
    lamda_s0(ic)=10^(log10(delta_ifv(ic))-(b*Mmin)+eps_global);
end


%% S2 MODEL
% SEISMOGENIC INDEX : FIXED
% B VALUE : FIXED (AS IN SHAPIRO ET AL. 2010)
b_shapiro=1.5;
eps_shapiro=0.25;
% Forecast using fixed b value and fixed seismogenic index, as given in
% Shapiro et al. 2010
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
         delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); 
    end
    lamda_s2(ic)=10^(log10(delta_ifv(ic))-(b_shapiro*Mmin)+eps_shapiro);
end


%% S3 MODEL
% SEISMOGENIC INDEX : FIXED
% B VALUE : FIXED (AS IN BACHMANN ET AL.)
b_corinne=1.0;
eps_corinne=(log10(rate_obs_update(23)))-(log10(ifv2(23)))+(b_corinne*Mcat)
% Forecast using fixed b value and fixed seismogenic index, as given in
% Shapiro et al. 2010
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
         delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); 
    end
    lamda_s3(ic)=10^(log10(delta_ifv(ic))-(b_corinne*Mmin)+eps_corinne);
end




% bprob=bv(ic);



%% S4 MODEL
% SEISMOGENIC INDEX : FREE
% B VALUE : FIXED (AS IN BACHMANN ET AL.)
b_corinne=1.0;
% Update seismogenic index every 6, 12, 18, 24, ....hours
eps_upd3(1)=0.5;
for ip=2:length(ifv2);
    eps_upd3(ip)=log10(rate_obs_update(ip))-log10(ifv2(ip))+b_corinne*Mcat; %SEISMOGENIC INDEX UPDATED AT 6, 12, 18 ... HOURS
end
% Forecast using fixed b value and fixed seismogenic index, as given in
% Shapiro et al. 2010
for ic=1:length(delta_ifv);
    if delta_ifv(ic)==0;
         delta_ifv(ic)=1;
    else
        delta_ifv(ic)=delta_ifv(ic); 
    end
    lamda_s4(ic)=10^(log10(delta_ifv(ic))-(b_corinne*Mmin)+eps_upd3(ic));
end




% %% PROBABILITY CALCULATION BASED ON S4 MODEL
% %% Probability that events with magnitude higher than mag would occur during fluid injection
% mag1=1;
% for ic=1:length(delta_ifv);
%     P1(ic)=exp(-delta_ifv(ic)*(10^(eps_upd3(ic)-b_corinne*mag1))); 
% end
% PP1=1-P1;
% %% Probability that events with magnitude higher than mag would occur during fluid injection
% mag2=2;
% for ic=1:length(delta_ifv);
%     P2(ic)=exp(-delta_ifv(ic)*(10^(eps_upd3(ic)-b_corinne*mag2))); 
% end
% PP2=1-P2;
% %% Probability that events with magnitude higher than mag would occur during fluid injection
% mag3=3;
% for ic=1:length(delta_ifv);
%     P3(ic)=exp(-delta_ifv(ic)*(10^(eps_upd3(ic)-b_corinne*mag3))); 
% end
% PP3=1-P3;
% %% Probability that events with magnitude higher than mag would occur during fluid injection
% mag4=4;
% for ic=1:length(delta_ifv);
%     P4(ic)=exp(-delta_ifv(ic)*(10^(eps_upd3(ic)-b_corinne*mag4))); 
% end
% PP4=1-P4;


% PROBABILITY CALCULATION BASED ON SR MODEL
% Probability that events with magnitude higher than mag would occur during fluid injection
mag1=1;
for ic=1:length(delta_ifv);
    P1(ic)=exp(-delta_ifv(ic)*(10^(eps_upd2(ic)-bv(ic)*mag1))); 
end
PP1=1-P1;
% Probability that events with magnitude higher than mag would occur during fluid injection
mag2=2;
for ic=1:length(delta_ifv);
    P2(ic)=exp(-delta_ifv(ic)*(10^(eps_upd2(ic)-bv(ic)*mag2))); 
end
PP2=1-P2;
% Probability that events with magnitude higher than mag would occur during fluid injection
mag3=3;
for ic=1:length(delta_ifv);
    P3(ic)=exp(-delta_ifv(ic)*(10^(eps_upd2(ic)-bv(ic)*mag3))); 
end
PP3=1-P3;
% Probability that events with magnitude higher than mag would occur during fluid injection
mag4=4;
for ic=1:length(delta_ifv);
    P4(ic)=exp(-delta_ifv(ic)*(10^(eps_upd2(ic)-bv(ic)*mag4))); 
end
PP4=1-P4;






% PLOTING
%
%Plot the entire catalog data
figure;
subplot(3,1,1);
plot(ctime,mag,'+');
title('The entire catalog data');
ylabel('Magnitude');
% Plot the first 6 days of the catalog
subplot(3,1,2);
plot(t15,mag15,'+');
title('15 days catalog');
ylabel('Magnitude');
% Plot events with magnitudes >= Mmin
subplot(3,1,3);
plot(t,MM,'+');
title('15 days catalog, M>=Mmin');
xlabel('Time (days)');
ylabel('Magnitude');  


mg=0.9:0.1:3;
for ii=1:length(mg);
     rt=MM>=mg(ii);
    ind=find(rt); 
    NN(ii)=length(ind);
    clear rt, ind
end



mg2=0:0.1:3;













figure; 
semilogy(mg,NN,'k.','MarkerSize',20); hold on;
Nnew=(10^fAValue).*(10.^(-b.*mg));
semilogy(mg,Nnew,'k','LineWidth',2);
ylabel('log(N)','FontSize',20);
xlabel('Magnitude','FontSize',20);
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','normal','FontName','Times');
grid on;


figure;
subplot(3,2,1);
plot(t15,mag15,'k.','MarkerSize',18);
ylabel('M_w','FontSize',18);
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times');
% grid on;

subplot(3,2,3);
plot(time,flow_rate,'k'); grid on;
ylabel('Flow Rate(l/min)','FontSize',18);
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times');
xlim([0 15]);

subplot(3,2,5);
plot(timet,(mm./10^7)./10000,'k','LineWidth',2);
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times');
set(get(gca,'XLabel'),'String','Time (days)','FontSize',18,'FontName','Times')
set(get(gca,'YLabel'),'String','M_0 (Nm/10^5)','FontSize',18,'FontName','Times');

subplot(1,2,2);
semilogy(mg,NN,'k.','MarkerSize',18); hold on;
Nnew=(10^fAValue).*(10.^(-b.*mg));
semilogy(mg,Nnew,'k','LineWidth',2);
ylabel('Cumulative Number','FontSize',18);
xlabel('M_w','FontSize',18);
grid on;
text(2.0,800,'a = 4.29','FontSize',18); hold on;
text(2.0,400,'b = 1.47','FontSize',18);
set(gca,'LineWidth',1,'FontSize',18,'FontWeight','normal','FontName','Times');




figure;
t1=fForecastStart:fDt:fForecastEnd;
subplot(2,1,1);
plot(t1,b_upd,'k'); grid on;
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','normal','FontName','Times');
ylim([1 2]);
ylabel('b value','FontSize',20);

subplot(2,1,2);
plot(t1,eps_upd2,'k'); grid on;
f=texlabel('epsilon');
ylabel(f,'FontSize',20);
xlabel('Time (days)','FontSize',20);
ylim([-1 1]);
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','normal','FontName','Times');







% catalog.time=ctime;
% catalog.t15=t15;
% catalog.mag=mag;
% catalog.mag15=mag15;
% 
% 
% savefile='catalog.mat';
% save(savefile,'catalog');

% % Ploting the flow rate, the injected fluid volume, and the observed seismicity rate
figure;
 t2=fForecastStart:fDt:fForecastEnd+fDt;
subplot(3,1,1);
plot(time,flow_rate); grid on;
ylabel('Flow Rate(l/min)','FontSize',16);
subplot(3,1,2);
plot(t2,ifv); grid on;
ylabel('Inj. Fluid Vol.(m^3)','FontSize',16);
subplot(3,1,3);
semilogy(t2,rate_obs,'ko'); grid on;
hold on;
semilogy(t2,rate_obs,'r'); 
title('(Mmin=0.9)','FontSize',16);
ylabel('Observed Rate','FontSize',16);
ylim([0.05 10^3]);
xlabel('Time (days)','FontSize',16);


% ftime=0.25:0.25:6;
% for ii=1:length(ftime);
%     xx=ftime(ii)-time;
%     xy=xx(xx>0);
%     tindic=find(xx==min(xy));
%     tind(ii)=min(tindic)
%     
% %     Frate(ii)=flow_rate(tind);
%     clear tindic;
% end
% 
% mean_flow_rate(1)=mean(flow_rate(1:tind(1)));
% for j=1:length(tind)-1;
%     mean_flow_rate(j+1)=mean(flow_rate(tind(j):tind(j+1)));
% end
%         
% 
% figure;
% subplot(2,1,1);
% plot(time,flow_rate); 
% xlim=([0 6]);
% subplot(2,1,2);
% plot(ftime,mean_flow_rate);
% 
% save mean_flow_rate.mat;

    
    




% 
% figure;
% subplot(2,1,1);
% plot(t,M,'+');
% ylabel('Magnitude','FontSize',14);  
% subplot(2,1,2);
% semilogy(t2,rate_obs,'ko'); grid on;
% hold on;
% semilogy(t2,rate_obs,'r'); 
% ylabel('Observed Rate','FontSize',14);
% ylim([0.05 10^3]);
% xlabel('Time (days)','FontSize',14);
% 
% 
% figure;
% t2=fForecastStart:fDt:fForecastEnd+fDt;
% subplot(2,1,1);
% plot(time,flow_rate); grid on;
% ylabel('Flow Rate(l/min)','FontSize',14);
% subplot(2,1,2);
% plot(t2,ifv); grid on;
% ylabel('Inj. Fluid Vol.(m^3)','FontSize',14);
% xlabel('Time (days)','FontSize',14);




% Ploting the seismogenic index
figure;
% subplot(2,1,1);
% plot(t2(2:end),delta_eps); grid on;
% ylabel('Seismogenic index');
% title('Within 6 hours bins');
% ylim([-1 3]);
% subplot(2,1,2);
plot(t2(2:end),eps_upd); grid on;
ylabel('Seismogenic index');
xlabel('Time (days)');
title('Updated every 6 hours');
ylim([-1 3]);


% % Ploting the seismogenic index
% figure;
% subplot(3,1,1);
% plot(tb,b_upd); grid on;
% ylabel('b value');
% title('Updated b value (for Basel sequence, S2 Model)');
% subplot(3,1,2);
% plot(t2(2:end),eps_upd); grid on;
% ylabel('Seismogenic index');
% title('Updated using fixed b = 1.47 (S1 Model)');
% ylim([-1 3]);
% subplot(3,1,3);
% plot(t2(2:end),eps_upd2); grid on;
% ylabel('Seismogenic index');
% xlabel('Time (days)');
% title('Updated using free b value (S2 Model)');
% ylim([-1 3]);
% % subplot(4,1,4);
% % plot(t2(2:end),eps_upd3); grid on;
% % ylabel('Seismogenic index');
% % xlabel('Time (days)');
% % title('Updated using fixed b = 1 (S5 Model)');
% % ylim([-1 3]);





% Ploting the rates within 6 hours bins
figure;
semilogy(t2,rate_obs,'k','LineWidth',8); grid on;
hold on;
semilogy(t2(2:end),lamda_s0,'Color',[0.4 0.4 0.4],'LineWidth',10);
hold on;
semilogy(t2(2:end),lamda_s1,'b','LineWidth',4);
hold on;
semilogy(t2(2:end),lamda_s2,'c','LineWidth',3);
% hold on;
% semilogy(t2(2:end),lamda_s3,'m','LineWidth',2);
hold on;
semilogy(t2(2:end),lamda_sr,'r','LineWidth',2);
legend('observed','S0','S1','S2','SR');
ylabel('Rates within 6 hourS');
xlabel('Time (days)');
xlim([0 10]);
set(gca,'LineWidth',1,'FontSize',24,'FontWeight','normal','FontName','Times');
set(get(gca,'XLabel'),'String','Time (days)','FontSize',24,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within 6 hours','FontSize',24,'FontName','Times');





% figure;
% plot(t2(2:end),P);
% 
% 
% % Ploting the rates within 6 hours bins
% figure;
% plot(t2,rate_obs,'k','LineWidth',4); grid on;
% hold on;
% % plot(t2(2:end),lamda_s0,'g','LineWidth',8);
% % hold on;
% % plot(t2(2:end),lamda_s1,'b','LineWidth',4);
% % hold on;
% % plot(t2(2:end),lamda_s2,'c','LineWidth',3);
% % hold on;
% % plot(t2(2:end),lamda_s3,'m','LineWidth',2);
% % hold on;
% plot(t2(2:end),lamda_sr,'r','LineWidth',2);
% legend('observed','SR');
% ylabel('Rates within 6 hourS');
% xlabel('Time (days)');
% xlim([0 10]);
















% Ploting the probability
figure;
plot(t2(2:end),PP1,'k'); hold on;
plot(t2(2:end),PP2,'r'); hold on;
plot(t2(2:end),PP3,'b'); hold on;
plot(t2(2:end),PP4,'g'); 
legend('M>1','M>2','M>3','M>4');
xlabel('Time (days)');
ylabel('Probability');
title('Probability of events within 6 hours, based on SR');
xlim([0 10]);


tdur = 5.75; % duration of fluid injection in days

%% Prepare forecast_file.mat as input into log likelihood calculation
% column 1: [minimum longitude of the spatial bin]
% column 2: [maximum longitude of the spatial bin]
% column 3: [minimum latitude]
% column 4: [max latitude]
% column 5: [0 (depth min)]
% column 6: [30 (depth max)]
% column 7: [min magnitude]
% column 8: [max magnitude]
% column 9: [rate forecast]
% TAKE DATA UNTIL THE FLUID INJECTION IS STOPPED

longmin=zeros(tdur/fDt,1);
longmax=zeros(tdur/fDt,1);
latmin=zeros(tdur/fDt,1);
latmax=zeros(tdur/fDt,1);
depth=zeros(tdur/fDt,1);
depthmax=zeros(tdur/fDt,1);
Mwmin=zeros(tdur/fDt,1);
Mwmax=zeros(tdur/fDt,1);
column10(:,1)=zeros(tdur/fDt,1);


lamda_s0=lamda_s0';
lamda_s1=lamda_s1';
lamda_s2=lamda_s2';
lamda_s3=lamda_s3';
lamda_sr=lamda_sr';
lamda_s00=lamda_s0(1:tdur/fDt);
lamda_s11=lamda_s1(1:tdur/fDt);
lamda_s22=lamda_s2(1:tdur/fDt);
lamda_s33=lamda_s3(1:tdur/fDt);
lamda_srr=lamda_sr(1:tdur/fDt);

longmin(:,1)=7.5904;
longmax(:,1)=7.6014;
latmin(:,1)=47.5786;
latmax(:,1)=47.5896;
depth(:,1)=0;
depthmax(:,1)=1;
Mwmin(:,1)=0.9;
Mwmax(:,1)=3;
column10(:,1)=1;

S0=[longmin longmax latmin latmax depth depthmax Mwmin Mwmax lamda_s00 column10];
S1=[longmin longmax latmin latmax depth depthmax Mwmin Mwmax lamda_s11 column10];
S2=[longmin longmax latmin latmax depth depthmax Mwmin Mwmax lamda_s22 column10];
S3=[longmin longmax latmin latmax depth depthmax Mwmin Mwmax lamda_s33 column10];
SR=[longmin longmax latmin latmax depth depthmax Mwmin Mwmax lamda_srr column10];

S0123R.S0=S0;
S0123R.S1=S1;
S0123R.S2=S2;
S0123R.S3=S3;
S0123R.SR=SR;

% save S0123R.mat;

lamda_SR=lamda_srr;
% save lamda_SR.mat;








































