
function [scpint3,scpint5,scpint7]=calc_int(b,data)





% load S0123R_15days_2.mat;  
% SR=S0123R_15days_2.SR; SRR=SR(1:60); SRR=SRR';
% data=SRR(17)+SRR(18)+SRR(19)+SRR(20);



%% Mmax=3.7
%from min mag to max mag

Mmin=0.9; Mmax=3.7; dm=0.1; r=0;
bgmax=[];
for m=Mmin:dm:Mmax;
    diff=Mmin-m;
    bgm=(data.*10.^(b*diff));
    bgmax=[bgmax bgm];
end

k=length(bgmax)+2;
mBack(:,3:k)=bgmax(:,:);

%mBack(:,1)=47.6;mBack(:,2)= 7.5; 
mBack(:,1)=1;mBack(:,2)= 4; %% arbitrarily selected location, corresponds to the r=0 distance
scpint3=[];
for nIntLevel=0:0.1:10;
     [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_banu(mBack,nIntLevel,Mmin,Mmax,dm,r);
    scpint3=[scpint3 scP];
    
end

%% Mmax=5
%from min mag to max mag
bgmax=[];
Mmin=0.9; Mmax=5.0;
for m=Mmin:0.1:Mmax;
    diff=0.9-m;
    bgm=(data.*10.^(b*diff));
    bgmax=[bgmax bgm];
end

k=length(bgmax)+2;
mBack(:,3:k)=bgmax(:,:);
mBack(:,1)=47.6;mBack(:,2)= 7.5; %% arbitrarily selected location, corresponds to the r=0 distance

scpint5=[];
for nIntLevel=0:0.1:10;
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_banu(mBack,nIntLevel,Mmin,Mmax,dm,r);  
    scpint5=[scpint5 scP];
    
end

%% Mmax=7
%from min mag to max mag
bgmax=[];
Mmin=0.9; Mmax=7.0;
for m=Mmin:0.1:Mmax;
    diff=0.9-m;
    bgm=(data.*10.^(b*diff));
    bgmax=[bgmax bgm];
end

k=length(bgmax)+2;
mBack(:,3:k)=bgmax(:,:);
mBack(:,1)=47.6;mBack(:,2)= 7.5; %% arbitrarily selected location, corresponds to the r=0 distance

scpint7=[];
for nIntLevel=0:0.1:10;
    [scgx,scgy,scP,logP] = makeSwissBackgroundHazGrid_banu(mBack,nIntLevel,Mmin,Mmax,dm,r);
    scpint7=[scpint7 scP];

end



