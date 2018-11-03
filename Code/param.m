function param(cond, type)
global p phen_type
if cond == 0
    p.growth = 0;
else p.growth = 0.0055;
end
%%extra constant or function
p.clpxp=1%um
%%parameters of cckap's function
p1 =  -3.491e-06;
p2 =   0.0007817;
p3 =    -0.04119 ;
p4 =      0.7378;
T=150;%period of Caulobacter
  t_d=rem(t,T); %return remainder after division t/T
p.cckap=p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4;

%% Synthesis and degradation rate constants [units --> 1/min]
p.ks_cpdr=
p.kd_cpdr=
