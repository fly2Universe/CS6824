function initial_cond = get_init(type)
%% INITIAL CONDITIONS
y0=zeros(501,1);
y0(501) = 0.013;           % cell size

%%pass parameter
param(0, eval(type));
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 150;%150;        
[t,y]=ode15s(@odes,[t0 tf],y0);

initial_cond.init = y(end,:);