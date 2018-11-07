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





end