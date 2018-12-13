clc
clear all
%%CpdR~p
 a=[0,20,40,60,80,100,120];
o=[31106.238,3300.255,1194.719,751.82,17808.832,14729.225,188039.817];
b=o/max(o);
values = spcrv([[a(1) a a(end)];[b(1) b b(end)]],3);
plot(values(1,:),values(2,:), 'k');
hold on;
%%CpdR
 a=[0,20,40,60,80,100,120];
 o=[19774.782,204118.296,13489.933,10158.74,10186.033,13034.447,20443.276];
b=o/max(o);
values = spcrv([[a(1) a a(end)];[b(1) b b(end)]],3);
plot(values(1,:),values(2,:), 'b');
hold on;

a=[0,20,40,60,80,100,120,140];
 o=[149.2456667,233.478,150.655,105.8717,89.15376667,97.50456667,82.6268,87.39566667];
b=o/max(o);
values = spcrv([[a(1) a a(end)];[b(1) b b(end)]],3);
plot(values(1,:),values(2,:), 'r');
legend('CpdR~p','CpdR','PdeA')
figure (2)
  a=[10,30,50,70,90,110,130];
  o=[21418.761,2987.841,816.82,1585.548,2981.669,11241.64,11459.083];
  b=o/max(o);
values = spcrv([[a(1) a a(end)];[b(1) b b(end)]],3);
plot(values(1,:),values(2,:), 'k');
legend('TacA')

