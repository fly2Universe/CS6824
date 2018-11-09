%ballfall.m
%simulate the bouncing of a free falling ball
clear;
t0=0;
tmax=100;
tspan=[t0,tmax];

y0=[10,0];
teout = [];
yeout = [];
ieout = [];

options=odeset('Events',@bounceEvents);
[t,y,te,ye,ie] = ode45(@ballode,tspan,y0,options);
nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];
   teout = [teout; te];          % Events at tstart are never reported.
   yeout = [yeout; ye];
   ieout = [ieout; ie];


%plot
figure;
plot(t,y(1));

%ballode.m

function dydt=ballode(~,y)
dydt=[y(2);-9.8];
end

%bounceEvents.m

function [value,isterminal,direction]= bounceEvents(~,y)
value = y(1);
isterminal = 1;
direction = -1;
end

