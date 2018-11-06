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
p.k1_pos=0.5 %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=0.2 %backward reaction: Complex1->ClpXP+CpdR
p.ks_cpdr=0.02;
p.kd_cpdr=0.01;
%%phosporylation
p.k2_pos=0.5;
p.k2_neg=0.8;
%% Diffusion parameters [units --> um^2/min]
p.D_complex1=10;
p.D_cpdr=100;
p.D_cpdrp=100;
%%
p.J1=0.5;
p.J2=0.5;
p.J3=0.5;
%%binding
p.kcpdr_b_f=0.1;
p.kcpdr_f_b=1;
p.kcpdr_b_f=0.1;
p.kcpdr_f_b=1;
