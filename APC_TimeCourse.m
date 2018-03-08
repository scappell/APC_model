%PROGRAM SIMULATIING APC/C ACTIVITY TIME-COURSES SHOWN IN FIG 4e in CAPPELL ET AL, NATURE 2018
%---T.Meyer 
clear all;close all;clc
%----Fixed parameters for APC/C inactivation model--------------------------------------------
AT = 500; %Estimated total concentration of APC/C (nM)
A = 3.7 ; %Max rate of ubiquitin elongation, binding plus elongation, assumes distributive, per minute
A0 = 10 ; %Rate of mono-ubiquitination of Emi1, overall much slower than polyubiquitination 
    %since only a fraction of non-ubiquitinated Emi1 is unbound, per minute
B = 2; %Rate of deubiquitination of Emi1 (gamma- term) is assumed to be fast compared to B and can 
    %therefore be incorporated into parameter B, per minute (assumed that B=B0)
K = 5 ; %Binding affinity of Emi1 to inhibitory site, nM
E = .0025 ;  %Slow APC-independent Emi1 degradation, 13 hours (800 min) turnover rate,
    %needed to prevent Emi1 levels to increase to infinity when APC is completely inhibited

%-----Parameters describing CDK2-mediated inhibition and Emi1 synthesis rate------
S = 10;  %Maximal Emi1 synthesis rate, ~20 nM/min, S(t/T2), E2F-driven increase in Emi1 mRNA levels 
T2=500;  %Time constant for E2F-stimulated increase in Emi1 mRNA
n=3;     %Inactivation of APC/C by CDK2 activity, cooperativity and time constant of CDK2 increase
T1=600;  %Time constant fro CDK2-mediated APC/C inactivation   

%--------------------------------------------------------------------------
%Steady state analysis of Emi1 synthesis rate versus Emi1 concentration
%Plot tests how increasing CDK2 gradually inhibits APC, cooperativity of 3
%is assumed; 
tspan = [0 900];  %Simulation over 900 minutes
y0 = [100 0 0 0 0 0];  %Start parameters, y(1)-y(6), Emi1 concentrations for species with 0-5 ubiquitins added

[t3,y] = ode45(@(t,y) APC1a(t,y,A,B,A0,S,K,E,AT,n,T1,T2), tspan, y0);
dy=y(2:end,1)-y(1:end-1,1);
dt=t3(2:end)-t3(1:end-1);
figure(1),hold on

%Calculates APC activity from the total concentration of non-ubiquitinated
%Emi1 and CDK2 activity
z1=1./(1+(t3/T1).^n) .*(AT-((K+AT+y(:,1))/2 - 0.5*sqrt((K+AT+y(:,1)).^2 - 4*y(:,1)*AT)));
Norm=max(z1);
plot(t3(t3>50)/60,z1(t3>50)/Norm,'k','linewidth',3)

[t3,y] = ode45(@(t,y) APC1b(t,y,A,B,A0,S,K,E,AT,T2), tspan, y0);
dy=y(2:end,1)-y(1:end-1,1);
dt=t3(2:end)-t3(1:end-1);

%Calculates APC activity from the total concentration of non-ubiquitinated Emi1
%without CDK2 activity
z2=(AT-((K+AT+y(:,1))/2 - 0.5*sqrt((K+AT+y(:,1)).^2 - 4*y(:,1)*AT))); 
plot(t3(t3>50)/60,z2(t3>50)/Norm,'b','linewidth',3),hold on

%Calculates APC activity if there is no Emi1 increase
z3=1./(1+(t3/T1).^n);
plot(t3(t3>50)/60,z3(t3>50),'g','linewidth',3),hold on

axis([0 15 0 1.05])
xlabel('Time (hrs)','fontsize',16);
ylabel('Relative APC/C activity','fontsize',16);
set(gca,'fontsize',16,'linewidth',2,'tickdir','out');
title('Kinetic Simulation','fontsize',16);
legend({'WT','Cyclin E absent','Emi1 absent'},'box','off','fontsize',16);