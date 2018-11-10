%main_temp_odes
clc
clear all
close all
global p
param(1); %cell can grow
% p.growth = 0.0055;
% p.clpxp=1;%um
% p.k1_pos=0.5; %forward reaction: ClpXP+CpdR ->Complex1
% p.k1_neg=0.2; %backward reaction: Complex1->ClpXP+CpdR
% p.ks_cpdr=2;%0.02;
% p.kd_cpdr=1;%0.01;
% %%phosporylation
% p.k2_pos=2;%0.5;
% p.k2_neg=8;%0.8;%based on ratio of p/unp
% p.D_complex1=10;
% p.D_cpdr=100;
% p.D_cpdrp=100;
% p.J1=0.5;
% p.J2=0.5;
% p.J3=0.1;
% %%binding
% p.kcpdr_b_f=0.1;
% p.kcpdr_f_b=1;
% p.kcpdrp_b_f=0.5;
% p.kcpdrp_f_b=0.5;


[T,Y] = ode15s('temp_odes',[0,800],[1,1,1,1,1]);
figure(100);
subplot(3,1,1);
plot(T,Y(:,1),'g');
xlabel('time')
ylabel('species')
legend('p1')
subplot(3,1,2);
plot(T,Y(:,5),'b');
xlabel('time')
ylabel('species')
legend('CpdRp_b')
subplot(3,1,3);
plot(T,Y(:,3),'r');
xlabel('time')
ylabel('species')
legend('CpdR_b')