%PROGRAM SIMULATIING STEADY STATE relationship between rate of Emi1 concentration change as a
%function of Emi1 concentration before and after partial CDK2 inhibition, FIG 3i in CAPPELL ET AL, NATURE 2018
%---T.Meyer 
clear all;close all;clc;
AT = 500; %Estimated total concentration of APC/C (nM)
A = 3.7 ; %Max rate of ubiquitin elongation, binding plus elongation, assumes distributive, per minute
A0 = 10 ; %Rate of mono-ubiquitination of Emi1, overall much slower than polyubiquitination 
    %since only a fraction of non-ubiquitinated Emi1 is unbound, per minute
B = 2; %Rate of deubiquitination of Emi1 (gamma- term) is assumed to be fast compared to B and can 
    %therefore be incorporated into parameter B, per minute (assumed that B=B0)
K = 5 ; %Binding affinity of Emi1 to inhibitory site, nM
E = .0025 ;  %Slow APC-independent Emi1 degradation, 13 hours (800 min) turnover rate,
    %needed to prevent Emi1 levels to increase to infinity when APC is completely inhibited

S=3.5;
Inh=1;

%-------------Generates three plots below-------------- 
%1. Stable steady state low Emi1 level when Emi1 starts from a low concentration
%2. Switch to high Emi1 level when E2F is on and APC is partially inhibited by CDK2
%3. Cell remains at high Emi1 level even after CDK2 is again inhibited
tspan = [0 900];
y0 = [0 50 0 0 0 0];

[t3,y] = ode45(@(t,y) APC1c(t,y,A,B,A0,S,K,E,AT,Inh), tspan, y0);
yL=y(end,1);
y0 = [1200 50 0 0 0 0];
[t3,y] = ode45(@(t,y) APC1c(t,y,A,B,A0,S,K,E,AT,Inh), tspan, y0);
yH=y(end,1);

Inh=1;
figure(1),hold on
for i=1:500
    E0=i*4;
    tspan = [0 900];
    y0 = [E0 0 0 0 0 0];
    [t3,y] = ode45(@(t,y) APC1c(t,y,A,B,A0,S,K,E,AT,Inh), tspan, y0);
    n0=find(t3>20,1,'first');
    C(i)=log10(y(n0,1));
    R(i)=(y(n0,1)-y(n0-10,1))/(t3(n0)-t3(n0-10));
end
plot(C,R,'k-','linewidth',2)
v=hline(0,'k');
scatter(log10(yL),0,100,'r','filled')
axis([1.6 3.5 -1.7 2.7 ])
ylabel('dEmi1/dt, synthesis/degradation rate (nM/min)','fontsize',16)
xlabel('Emi1 concentration (log10, nM)','fontsize',16)
title('Stable low level of Emi1 when Emi1 levels start from low','fontsize',16)
set(gca,'fontsize',16,'linewidth',2,'tickdir','out');

Inh=0.85;
figure(2),hold on
for i=1:500
    E0=i*4;
    tspan = [0 900];
    y0 = [E0 0 0 0 0 0];
    [t3,y] = ode45(@(t,y) APC1c(t,y,A,B,A0,S,K,E,AT,Inh), tspan, y0);
    n0=find(t3>20,1,'first');
    C(i)=log10(y(n0,1));
    R(i)=(y(n0,1)-y(n0-10,1))/(t3(n0)-t3(n0-10));
end
plot(C,R,'k-','linewidth',2)
v=hline(0,'k');
scatter(log10(yH),0,100,'r','filled')
axis([1.6 3.5 -1.7 2.7 ])
ylabel('dEmi1/dt, synthesis/degradation rate (nM/min)','fontsize',16)
xlabel('Emi1 concentration (log10, nM)','fontsize',16)
title('Switch to high Emi1 when CDK2 partially inhibits APC','fontsize',16)
set(gca,'fontsize',16,'linewidth',2,'tickdir','out');

Inh=1;
figure(3),hold on

for i=1:500
    E0=i*4;
    tspan = [0 900];
    y0 = [E0 0 0 0 0 0];
    [t3,y] = ode45(@(t,y) APC1c(t,y,A,B,A0,S,K,E,AT,Inh), tspan, y0);
    n0=find(t3>20,1,'first');
    C(i)=log10(y(n0,1));
    R(i)=(y(n0,1)-y(n0-10,1))/(t3(n0)-t3(n0-10));
end
plot(C,R,'k-','linewidth',2)
v=hline(0,'k');
scatter(log10(yH),0,100,'r','filled')
axis([1.6 3.5 -1.7 2.7 ])
ylabel('dEmi1/dt, synthesis/degradation rate (nM/min)','fontsize',16)
xlabel('Emi1 concentration (log10, nM)','fontsize',16)
title('Cell remains in high Emi1 state even when CDK2 activity drops again','fontsize',16)
set(gca,'fontsize',16,'linewidth',2,'tickdir','out');

%% Function to draw horizontal line
function hhh=hline(y,in1,in2)
% function h=hline(y, linetype, label)
% 
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    x=get(gca,'xlim');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
        else
            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
        end
    end

    if g==0
    hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end
end