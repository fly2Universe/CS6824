%main_temp_odes
clc
clear all
close all
global p
param(1); %cell can grow


[T,Y] = ode15s('temp_odes',[0,500],[2,1,4,1,1]);
figure(100);
subplot(3,1,1);
plot(T,Y(:,1),'g');
xlabel('time')
ylabel('species')
legend('complex1')
subplot(3,1,2);
plot(T,Y(:,5),'b');
xlabel('time')
ylabel('species')
legend('complex2')
subplot(3,1,3);
plot(T,Y(:,6),'r');
xlabel('time')
ylabel('species')
legend('RcdA')