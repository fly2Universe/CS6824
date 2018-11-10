%temperal_odes.m

function dydt=temp_odes(t,y)
global p;
%%parameters of cckap's function
p1 =  -3.491e-06;
p2 =   0.0007817;
p3 =    -0.04119 ;
p4 =      0.7378;
T=120;%period of Caulobacter
t_d=rem(t,T); %return remainder after division t/T
p.cckap=p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4;

%dydt=zeros(601,1);
%P1=y(1); CpdR_f=y(2); CpdR_b=y(3); CpdR~p_f=y(4); CpdR~p_b=y(5)
%%Complex1 (ClpXP:CpdR)
  dydt(1)=p.k1_pos*p.clpxp*y(3)-p.k1_neg*y(1);
%%CpdR_f 
 % dydt(2)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(1)*y(2)/(y(1)+p.J2)...
  %    +p.k2_pos*y(4)-p.k2_neg*p.cckap*y(2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(3)-p.kcpdr_f_b*y(2)
 dydt(2)=p.ks_cpdr*p.J1/(y(1)+p.J1)-p.kd_cpdr*y(1)*y(2)/(y(1)+p.J2)...
    +p.k2_pos*y(4)-p.k2_neg*p.cckap*y(2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(3)-p.kcpdr_f_b*y(2)
%%CpdR_b
  dydt(3)=-p.kd_cpdr*y(1)*y(3)/(y(1)+p.J2)+p.kcpdr_f_b*y(2)-p.kcpdr_b_f*y(3)...
      +p.k2_pos*y(5)-p.k2_neg*p.cckap*y(3)/(p.cckap+p.J3)+p.k1_neg*y(1)-p.k1_pos*y(3)*p.clpxp;

%%CpdR~P_f
  dydt(4)=p.k2_neg*p.cckap*y(2)/(p.cckap+p.J3)-p.k2_pos*y(4)+p.kcpdrp_b_f*y(5)-p.kcpdrp_f_b*y(4)

%%CpdR~P_b
dydt(5)=p.k2_neg*p.cckap*y(3)/(p.cckap+p.J3)-p.k2_pos*y(5)+p.kcpdrp_f_b*y(4)-p.kcpdrp_b_f*y(5);

%%%%%%%%%%%%%%%%%%%%%%
%dydt(6)=p.growth*y(6);
%%%%
dydt=dydt';
%%%%%%%%%%%add model 2