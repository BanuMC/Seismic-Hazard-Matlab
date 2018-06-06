clc
clear all
close all

load com.mat;
% data=com(22);

datax=sum(com(1:24));


b2=1.0;
[scp3,scp5,scp7]=calc_int(b2,datax)

%% PLOTING
Int=0:0.1:10;
figure;
semilogy(Int,scp3,'LineWidth',10,'Color',[.9 .9 .9]); hold on;
semilogy(Int,scp5,'LineWidth',6,'Color',[.7 .7 .7]); hold on;
semilogy(Int,scp7,'LineWidth',2,'Color',[.1 .1 .1]); hold on;
ylabel('Probability of exceeding EMS intensity','FontSize',20);
xlabel('EMS Intensity','FontSize',20);
legend('Mmax=3.7','Mmax=5','Mmax=7');
ylim([10^-4 10^0]);
xlim([2 10]);
grid on
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','EMS Intensity','FontSize',24,'FontName','Times')
set(get(gca,'YLabel'),'String','Probability of exceeding EMS intensity','FontSize',24,'FontName','Times')



% %% PLOTING
% Int=0:0.1:10;
% figure;
% semilogy(Int,scp5,'LineWidth',2,'Color',[.1 .1 .1]); hold on;
% ylabel('Probability of exceeding EMS intensity','FontSize',20);
% xlabel('EMS Intensity','FontSize',20);
% ylim([10^-4 10^0]);
% xlim([2 10]);
% grid on
% set(gca,'LineWidth',2,'FontSize',20,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','EMS Intensity','FontSize',20,'FontName','Times')
% set(get(gca,'YLabel'),'String','Probability of exceeding EMS intensity','FontSize',20,'FontName','Times')






% data1=com(2);
% data2=com(8);
% data3=com(12);
% data4=com(16);
% data5=com(20);
% data6=com(24);
% data7=com(36); %day 9
% data8=com(48); % day12
% data9=com(58); % day15
% 
% 
% 
% 
% Mmin=0.9; 
% Mmax=5;
% 
% [haz1]=calc_int_M(b,data1,Mmin,Mmax);
% [haz2]=calc_int_M(b,data2,Mmin,Mmax);
% [haz3]=calc_int_M(b,data3,Mmin,Mmax);
% [haz4]=calc_int_M(b,data4,Mmin,Mmax);
% [haz5]=calc_int_M(b,data5,Mmin,Mmax);
% [haz6]=calc_int_M(b,data6,Mmin,Mmax);
% [haz7]=calc_int_M(b,data7,Mmin,Mmax);
% [haz8]=calc_int_M(b,data8,Mmin,Mmax);
% [haz9]=calc_int_M(b,data9,Mmin,Mmax);
% 
% figure;
% semilogy(Int,haz1,'k','LineWidth',2); hold on;
% semilogy(Int,haz2,'LineWidth',5,'Color',[.8 .8 .8]); hold on;
% semilogy(Int,haz3,'LineWidth',9,'Color',[.5 .5 .5]); hold on;
% semilogy(Int,haz4,'LineWidth',4,'Color',[.3 .3 .3]); hold on;
% semilogy(Int,haz5,'k--','LineWidth',2); hold on;
% semilogy(Int,haz6,'k-.','Color',[.6 .6 .6],'LineWidth',2); hold on;
% semilogy(Int,haz7,'k-.','LineWidth',2); hold on;
% semilogy(Int,haz8,'-.','LineWidth',2,'Color',[.4 .4 .4]); hold on;
% semilogy(Int,haz9,'-.','LineWidth',2,'Color',[.8 .8 .8]); hold on;
% % % 
% 
% 
% ylabel('Probability of exceeding EMS intensity','FontSize',20);
% xlabel('EMS Intensity','FontSize',20);
% legend('Day1','Day2','Day3','Day4','Day5','Day6','Day9','Day12','Day15');
% % ylim([10^-3 10^0]);
% % xlim([2 7]);
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','EMS Intensity','FontSize',24,'FontName','Times')
% set(get(gca,'YLabel'),'String','Probability of exceeding EMS intensity','FontSize',24,'FontName','Times')
% 
% 
% 
% 
% 















