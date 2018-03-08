function dydt = APC1a(t,y,A,B,A0,S,K,E,AT,n,T1,T2)
%dydt from 1 to 6 corresponds to the concentration change of the species of Emi1 without
%conjugated ubiquitin (y(1) corresponds to E0), to the species with 1 to 5 
%conjugated ubiquitins (y(2) to y(6) corresponds to E1 to E5); 
%It is assumed that the species with 5 conjugated ubiquitins is rapidly degraded
%The function AE calculates the concentration of APC inhibited y Emi1
%1/(1+(t/T1)^n) describes the time course of inhibition of APC by CDK2
%S*(1-exp(-t/T2)) describes the E2F-regulated increase in Emi1 mRNA
%It is assumed that Emi1 bound to the inhibitory sites on APC cannot be ubiquitinated and that
%ubiquitinated Emi1 cannot inhibit APC

dydt = [0 0 0 0 0 0]';
dydt(1) = S*(1-exp(-t/T2))  +  B* y(2) - A0/(1+(t/T1)^n)/AT*(AT-AE(y(1),AT,K)) *(y(1)-AE(y(1),AT,K))- E * y(1);  % y(1) total Emi1 without 
dydt(2) = -  B* (y(2)-y(3)) + A0/(1+(t/T1)^n)/AT*(AT-AE(y(1),AT,K)) *(y(1)-AE(y(1),AT,K)) - A/(1+(t/T1)^n)/AT*(AT-AE(y(1),AT,K)) * y(2);
dydt(3) = -  B* (y(3) -y(4)) + A/(1+(t/T1)^n)/AT* (AT-AE(y(1),AT,K)) * (y(2) -y(3));
dydt(4) = -  B* (y(4) -y(5)) + A/(1+(t/T1)^n)/AT*  (AT-AE(y(1),AT,K)) * (y(3) -y(4));
dydt(5) = -  B* (y(5) -y(6)) + A/(1+(t/T1)^n)/AT*  (AT-AE(y(1),AT,K)) * (y(4) -y(5));
dydt(6) = -  B * y(6) + A/(1+(t/T1)^n)/AT* (AT-AE(y(1),AT,K)) * (y(5)-y(6));

function y=AE(y1,AT,K)  %Function calculating the concentration of inhibited APC-Emi1 complexes
y=0.5*((K+AT+y1) - sqrt((K+AT+y1)^2 - 4*y1*AT));