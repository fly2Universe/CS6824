%odes.m


function dydt=odes(t,y)
global p;
%%parameters of cckap's function
% p1 =  -3.491e-06;
% p2 =   0.0007817;
% p3 =    -0.04119 ;
% p4 =      0.7378;
T=120;%period of Caulobacter
t_d=rem(t,T); %return remainder after division t/T
%p.cckap=p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4;
p1 =   1.315e-06 ;
       p2 =  -0.0002346 ;
       p3 =    0.001233;
       p4 =      0.9975;
       p.divkp_free=p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4;
   %%%%DivK~P_oldpole
    d1 =  -8.515e-06;
       d2 =   0.0009505;
       d3 =    -0.01231 ;
       d4 =      0.1302 ;
     p.divkp_oldpole= (d1*t_d.^3 + d2*t_d.^2 + d3*t_d + d4).*(t_d<=98)+0.*(t_d>98);

dydt=zeros(601,1);

%%Complex1 (ClpXP:CpdR)
dydt(1)=p.k1_pos*p.clpxp*y(1+200)/(y(201)+p.km1)-p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(2)-y(1));
for v1=2:99
  dydt(v1)=p.k1_pos*p.clpxp*y(v1+200)/(y(v1+200)+p.km1)-p.k1_neg*y(v1)+(p.D_complex1/y(601)^2)*(y(v1-1)-2*y(v1)+y(v1+1));
end
dydt(100)=p.k1_pos*p.clpxp*y(100+200)/(y(300)+p.km1)-p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(99)-y(100));

%%CpdR_f 
% dydt(101)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
%       +p.k2_pos*y(101+200)-p.k2_neg*p.cckap*y(101)/(p.cckap+p.J3)+p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)...
%       +(p.D_cpdr/y(601)^2)*(y(102)-y(101));
% for v2=102:199
%   dydt(v2)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
%       +p.k2_pos*y(v2+200)-p.k2_neg*p.cckap*y(v2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2)...
%       +(p.D_cpdr/y(601)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
% end
% dydt(200)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(200-100)*y(v2)/(y(200-100)+p.J2)...
%       +p.k2_pos*y(200+200)-p.k2_neg*p.cckap*y(200)/(p.cckap+p.J3)+p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)...
%       +(p.D_cpdr/y(601)^2)*(y(199)-y(200));
%%%%%%%%%%modified 11/10
% dydt(101)=p.ks_cpdr*p.J1/(y(101-100)+p.J1)-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
%       +p.k2_pos*y(101+200)-p.k2_neg*p.cckap*y(101)/(p.cckap+p.J3)+p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)...
%       +(p.D_cpdr/y(601)^2)*(y(102)-y(101));
% for v2=102:199
%   dydt(v2)=p.ks_cpdr*p.J1/(y(v2-100)+p.J1)-p.kd_cpdr*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
%       +p.k2_pos*y(v2+200)-p.k2_neg*p.cckap*y(v2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2)...
%       +(p.D_cpdr/y(601)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
% end
% dydt(200)=p.ks_cpdr*p.J1/(y(200-100)+p.J1)-p.kd_cpdr*y(200-100)*y(v2)/(y(200-100)+p.J2)...
%       +p.k2_pos*y(200+200)-p.k2_neg*p.cckap*y(200)/(p.cckap+p.J3)+p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)...
%       +(p.D_cpdr/y(601)^2)*(y(199)-y(200));
%%%%replace ccka with divk classified
dydt(101)=p.ks_cpdr*p.J1/(y(101-100)+p.J1)-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
      +p.k2_pos*y(101+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(101)+p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)...
      +(p.D_cpdr/y(601)^2)*(y(102)-y(101));
for v2=102:199
  dydt(v2)=p.ks_cpdr*p.J1/(y(v2-100)+p.J1)-p.kd_cpdr*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
      +p.k2_pos*y(v2+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(v2)+p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2)...
      +(p.D_cpdr/y(601)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
end
dydt(200)=p.ks_cpdr*p.J1/(y(200-100)+p.J1)-p.kd_cpdr*y(200-100)*y(v2)/(y(200-100)+p.J2)...
      +p.k2_pos*y(200+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(200)+p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)...
      +(p.D_cpdr/y(601)^2)*(y(199)-y(200));
%%CpdR_b
% for v3=201:300
%   dydt(v3)=-p.kd_cpdr*y(v3-200)*y(v3)/(y(v3-200)+p.J2)+p.kcpdr_f_b*y(v3+300)*y(v3-100)-p.kcpdr_b_f*y(v3)...
%       +p.k2_pos*y(v3+200)-p.k2_neg*p.cckap*y(v3)/(p.cckap+p.J3)+p.k1_neg*y(v3-200)-p.k1_pos*y(v3)*p.clpxp;
% end
for v3=201:300
  dydt(v3)=-p.kd_cpdr*y(v3-200)*y(v3)/(y(v3-200)+p.J2)+p.kcpdr_f_b*y(v3+300)*y(v3-100)-p.kcpdr_b_f*y(v3)...+p.k2_pos*y(v3+200)*p.divkp_oldpole/(p.divkp_oldpole+p.J4)-p.k2_neg*y(v3)
      +p.k1_neg*y(v3-200)-p.k1_pos*y(v3)*p.clpxp;
end

%%CpdR~P_f
% dydt(301)=p.k2_neg*p.cckap*y(101)/(p.cckap+p.J3)-p.k2_pos*y(301)+p.kcpdrp_b_f*y(401)-p.kcpdrp_f_b*y(301)...
% +p.D_cpdrp*(y(302)-y(301))/(y(601)^2);
% 
% for v4=302:399
%   dydt(v4)=p.k2_neg*p.cckap*y(v4-200)/(p.cckap+p.J3)-p.k2_pos*y(v4)+p.kcpdrp_b_f*y(v4+100)-p.kcpdrp_f_b*y(v4)...
% +p.D_cpdrp*(y(v4+1)-2*y(v4)+y(v4-1))/(y(601)^2);
% end
% dydt(400)=p.k2_neg*p.cckap*y(200)/(p.cckap+p.J3)-p.k2_pos*y(400)+p.kcpdrp_b_f*y(500)-p.kcpdrp_f_b*y(400)...
% +p.D_cpdrp*(y(399)-y(400))/(y(601)^2);
%%%%%replace ccka with divk
dydt(301)=p.k2_neg*y(101)-p.k2_pos*y(301)*p.divkp_free/(p.divkp_free+p.J3)+p.kcpdrp_b_f*y(401)-p.kcpdrp_f_b*y(301)...
+p.D_cpdrp*(y(302)-y(301))/(y(601)^2);

for v4=302:399
  dydt(v4)=p.k2_neg*y(v4-200)-p.k2_pos*y(v4)*p.divkp_free/(p.divkp_free+p.J3)+p.kcpdrp_b_f*y(v4+100)-p.kcpdrp_f_b*y(v4)...
+p.D_cpdrp*(y(v4+1)-2*y(v4)+y(v4-1))/(y(601)^2);
end
dydt(400)=p.k2_neg*y(200)-p.k2_pos*y(400)*p.divkp_free/(p.divkp_free+p.J3)+p.kcpdrp_b_f*y(500)-p.kcpdrp_f_b*y(400)...
+p.D_cpdrp*(y(399)-y(400))/(y(601)^2);
%%CpdR~P_b
% for v5=401:500
% dydt(v5)=p.k2_neg*p.cckap*y(v5-200)/(p.cckap+p.J3)-p.k2_pos*y(v5)+p.kcpdrp_f_b*y(v5-100)-p.kcpdrp_b_f*y(v5);
% end
%%%%divk
for v5=401:500
dydt(v5)=...p.k2_neg*y(v5-200)-p.k2_pos*y(v5)*p.divkp_oldpole/(p.divkp_oldpole+p.J4)+
    p.kcpdrp_f_b*y(v5-100)-p.kcpdrp_b_f*y(v5);
end
%%CpdR sticky
dydt(501:600)=0;

%%%%%%%%%%%%%%%%%%%%%%
dydt(601)=p.growth*y(601);
%%%%

%%%%%%%%%%%add model 2
