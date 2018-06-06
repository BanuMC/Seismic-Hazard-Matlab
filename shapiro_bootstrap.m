clc
clear all
% close all

load rate_obs15d.mat;


%% INPUT
fForecastEnd = 15; % End of the forecast window
Mcat=0.9;


cat=load (['catalog/','basel_spimed_25708.mat']);
[MM,t] =  cat2mag(cat,fForecastEnd,Mcat);




%% BOOTSTRAPPING THE DATA
% generate random number vector using a uniform distribution
% with values from 1 to length(x)
adum=1;
bdum=length(MM);

itmax=100;

for it=1:itmax;
    clear rdum xnew
    rdum = round(adum + (bdum-adum).*rand(length(MM),1)); % random number vector
    for samp=1:bdum
        MMnew(samp)=MM(rdum(samp));
        tnew(samp)=t(rdum(samp));
    end
    %calculate the model from  bootstrap data
    [lamda(it).sr,rate_obs] = model_SR(MMnew,tnew);
    [S_PI_upd(it).srpi] = model_SR_pi(rate_obs);
    lamda_all_sr(it).sr=[lamda(it).sr(1:23) S_PI_upd(it).srpi];
end


for ik=1:60;
    for it=1:itmax;
        X(ik,it)=lamda_all_sr(it).sr(ik);
    end
end

time=0.01:0.25:15;

tt=zeros(length(time),itmax); 
for ii=1:length(time);
    tt(ii,:)=time(ii);
end



for j=1:length(time);
    sigma(j)=std(X(j,:));
    avr(j)=mean(X(j,:));
end


load S0123R_15days_2.mat;  
SR=S0123R_15days_2.SR; SRR_model=SR(1:60); 
time=0.01:0.25:15;

figure;
% subplot(2,1,2);
plot(time,rate_obs15d(2:end),'r*','MarkerSize',10); hold on;
plot(time,SRR_model,'bo','MarkerSize',10);hold on;
errorbar(time,avr,sigma,'ks','MarkerSize',10);
set(gca,'LineWidth',1,'FontSize',24,'FontWeight','normal','FontName','Times');
set(get(gca,'XLabel'),'String','Time (days)','FontSize',24,'FontName','Times')
set(get(gca,'YLabel'),'String','Rates within six hours','FontSize',24,'FontName','Times');
axis([0 15 0 100]);
legend('Data','Model','Model Mean');










