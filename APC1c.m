function dydt = APC1c(t,y,A,B,A0,S,K,E,AT,Inh)

dydt = [0 0 0 0 0 0]';
dydt(1) = S  +  B* y(2) - Inh*A0/AT*(AT-AE(y(1),AT,K)) *(y(1)-AE(y(1),AT,K))- E * y(1);  % y(1) total Emi1 without 
dydt(2) = -  B* (y(2)-y(3)) + Inh*A0/AT*(AT-AE(y(1),AT,K)) *(y(1)-AE(y(1),AT,K)) - Inh*A/AT*(AT-AE(y(1),AT,K)) * y(2);
dydt(3) = -  B* (y(3) -y(4)) + Inh*A/AT* (AT-AE(y(1),AT,K)) * (y(2) -y(3));
dydt(4) = -  B* (y(4) -y(5)) + Inh*A/AT*  (AT-AE(y(1),AT,K)) * (y(3) -y(4));
dydt(5) = -  B* (y(5) -y(6)) + Inh*A/AT*  (AT-AE(y(1),AT,K)) * (y(4) -y(5));
dydt(6) = -  B * y(6) + Inh*A/AT* (AT-AE(y(1),AT,K)) * (y(5)-y(6));

function y=AE(y1,AT,K)
y=0.5*((K+AT+y1) - sqrt((K+AT+y1)^2 - 4*y1*AT));