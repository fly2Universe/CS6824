function initial_cond = get_init()
%% INITIAL CONDITIONS
y0=zeros(1101,1);
%y0(501:700)=1;
y0(1101) = 0.013;           % cell size

%%pass parameter
param(0);
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 300;%150;        
[t,y]=ode15s(@odes,[t0 tf],y0);

  figure(1)
  plot (t,y(:,1001))
  hold on;
 figure(2)
% plot(t,y(:,1))
 plot(t,y(:,701)+y(:,801))
 hold on;
  figure(3)
  plot(t,y(:,601))

initial_cond.init = y(end,:);