function initial_cond = get_init()
%% INITIAL CONDITIONS
y0=zeros(601,1);
y0(501) = 0.013;           % cell size

%%pass parameter
param(0);
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 200;%150;        
[t,y]=ode45(@odes,[t0 tf],y0);
plot (t,y(:,20))
%plot(t,y(:,501))
initial_cond.init = y(end,:);