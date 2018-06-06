function [b1, b, tb] =  calc_bvalue(Mcat,tmax,fDt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cat_name : catalog name
%%% Mmin : completeness magnitude (minimum magnitude)
%%% tmax : (fForecastEnd) End of the forecast window in days (for example 5.5 days)
%%% fDt : Time step (for examle 0.25 days corresponding to 6 hours time bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% b1 :  b value for the whole sequence
%%% b: b value vector for the updated catalog for every 6, 12, 18 hours

%% Load the catalog
cat=load (['catalog/','basel_spimed_25708.mat']);
aa=zeros(3560,10);
aa(:,(1:9))=cat.a;

%% start of injection
date_start=datenum(2006,12,2,18,00,00);
%% time in days
time = (datenum(floor(aa(:,3)),aa(:,4),aa(:,5),aa(:,8),aa(:,9),aa(:,10))-date_start); 
%% Magnitudes
M=aa(:,6);

%% Take 5.5 days only
t15=time(time<=tmax);
mag15=M(1:length(t15));
%% Events with magnitudes greated than or equal to the minumum magnitude
ind=find(mag15>=Mcat);
%% Time vector corresponding to events of M>=Mmin
for i=1:length(ind);
    t(i)=t15(ind(i));
    mag(i)=mag15(ind(i));
end
t=t';
mag=mag';

%% b value for the whole sequence of 5.5 days
[mean_m1, b1, sig1, av2] =  bmemag(mag);

%% Compute updated b value for every 6, 12, 18 24, ...  hours (0.25, 0.5, 0.75,... days)

for ij=1:(tmax/fDt);
%      indic=find(t(t>fDt & t<fDt*(ij+1)));
%      indic=find(t(t>fDt*2 & t<fDt*(ij+2)));
     indic=find(t(t>fDt*4 & t<fDt*(ij+4)));
     if length(indic)>0;
        for ik=1:length(indic);
            mm(ik)=mag(indic(ik));
        end
     else
         mm=0;
     end
    [mean_m(ij), b(ij), sig(ij), av(ij)] =  bmemag(mm);
    clear indic, mm;
end

tb=(fDt:fDt:tmax);
% 
% %% Ploting
% figure;
% plot(tb,b,'o');
% xlabel('Time (days)');
% ylabel('b value');




