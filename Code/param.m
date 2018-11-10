function param(cond)
global p
if cond == 0
    p.growth = 0;
else p.growth = 0.0055;
end
%%extra constant or function
p.clpxp=1;%um


%% Synthesis and degradation rate constants [units --> 1/min]
p.k1_pos=0.5; %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=0.2; %backward reaction: Complex1->ClpXP+CpdR
p.ks_cpdr=2;%0.02;
p.kd_cpdr=1;%0.01;
%%phosporylation
p.k2_pos=2;%0.5;
p.k2_neg=8;%0.8;%based on ratio of p/unp
%% Diffusion parameters [units --> um^2/min]
p.D_complex1=10;
p.D_cpdr=100;
p.D_cpdrp=100;
%%
p.J1=0.5;
p.J2=0.5;
p.J3=0.1%%0.5;
%%binding
p.kcpdr_b_f=0.1;
p.kcpdr_f_b=1;
p.kcpdrp_b_f=0.1;
p.kcpdrp_f_b=1;
