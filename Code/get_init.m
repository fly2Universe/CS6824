function initial_cond = get_init()
%% INITIAL CONDITIONS
y0=zeros(601,1);
y0(601) = 0.013;           % cell size

%%pass parameter
param(0);
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 300;%150;        
[t,y]=ode15s(@odes,[t0 tf],y0);
plot (t,y(:,1))
%plot(t,y(:,501))
initial_cond.init = y(end,:);