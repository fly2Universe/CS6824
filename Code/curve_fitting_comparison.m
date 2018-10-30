clc
clear all
%fit curve
 x = [10,20,40,60,80,100,120];
 o = [0.365958, 0.247904, 0.125544, 0.295375, 0.632501, 0.99639, 0.998939];
 c = o/max(o);
%cftool
%x as x-axis, c as y-axis, polynomial, degree 3
%%%%%%% %result% %%%%%%%%%%%%
% f(x) = p1*x^3 + p2*x^2 + p3*x + p4
% Coefficients (with 95% confidence bounds):
%        p1 =  -3.491e-06  (-5.432e-06, -1.55e-06)
%        p2 =   0.0007817  (0.0003982, 0.001165)
%        p3 =    -0.04119  (-0.06252, -0.01985)
%        p4 =      0.7378  (0.4382, 1.038)
% 
% Goodness of fit:
%   SSE: 0.008382
%   R-square: 0.9892
%   Adjusted R-square: 0.9784
%   RMSE: 0.05286

%%%%%%%%compare experimental data with fitting curve%%%%%%%%%%%
x=[x,x+150,x+300,x+450];
c=[c,c];
c=[c,c];
scatter(x,c,'b')
hold on;
T=150;%period
t1=0:1:600;
  t_d=rem(t1,T); %return remainder after division t/T
p1 =  -3.491e-06;
p2 =   0.0007817;
p3 =    -0.04119 ;
p4 =      0.7378;
 plot(t1,p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4)
 xlabel('time')
ylabel('CckA~P')
legend('experiment','simulation')

