function param(cond)
global p
if cond == 0
    p.growth = 0;
else p.growth = 0.0055;
end
%%extra constant or function
p.clpxp=0.5;%um


%% Synthesis and degradation rate constants [units --> 1/min]
p.k1_pos=2; %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=15; %backward reaction: Complex1->ClpXP+CpdR
p.ks_cpdr=2;%0.02;
p.kd_cpdr=1;%0.01;
%%phosporylation
p.k2_pos=2;%0.5;
p.k2_neg=15;%0.8;%based on ratio of p/unp
%%
p.k3_pos=20;%forward reaction: complex1+RcdA ->Complex1
p.k3_neg=0.2;%backward reaction: complex1+RcdA
p.ks_rcda=2.5;%0.023;
p.kd_rcda=2;%0.017;%from Tyson lab
%% Diffusion parameters [units --> um^2/min]
p.D_complex1=1;
p.D_cpdr=100;
p.D_cpdrp=100;
p.D_complex2=1;
p.D_rcda=100;
%%
p.J1=10;
p.J2=0.5;
p.J3=0.5;%%0.5;
p.J4=0.4;
p.km1=30;
p.J5=0.5;
p.J6=2;
%%binding
p.kcpdr_b_f=0.1;
p.kcpdr_f_b=15;
p.kcpdrp_b_f=0.5;
p.kcpdrp_f_b=1;
