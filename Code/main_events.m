%main_events.m

function sw_to_early_pd = main_events(initial_cond)
param(1); %cell can grow

%% INITIAL CONDITIONS
y0 = initial_cond.init;
y0 = y0.';
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 120;

options = odeset('Events',@event,'RelTol',1e-4,'AbsTol',1e-6);
tout=t0;
yout=y0.';
teout = [];
yeout = [];
ieout = [];

while t0<tf
    
    [t,y,te,ye,ie]=ode15s(@odes,[t0 tf],y0,options);
    
    nt=length(t);
    
    tout=[tout;t(2:nt)];
    yout=[yout;y(2:nt,:)];
    teout = [teout;te];
    yeout = [yeout;ye];
    ieout = [ieout;ie];
    
    y0 = y(nt,:);

    if isscalar(ie) == 0
        ie = 0;
    end
    
     t0=t(nt);
    %enforced localization of specific species
     
    if t0>=tf
        break;
    end
    
    sw_to_early_pd.endpoints = yout(end,:);

end