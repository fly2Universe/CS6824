%odes.m


function dydt=odes(t,y)
global p;

T=120;%period of Caulobacter
t_d=rem(t,T); %return remainder after division t/T

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

dydt=zeros(501,1);

%%%%%%%%%%Complex1 (ClpXP:CpdR)
dydt(1)=p.k1_pos*p.clpxp*y(1+200)/(y(201)+p.km1)-p.k1_neg*y(1)+(p.D_complex1/y(501)^2)*(y(2)-y(1));
for v1=2:99
  dydt(v1)=p.k1_pos*p.clpxp*y(v1+200)/(y(v1+200)+p.km1)-p.k1_neg*y(v1)+(p.D_complex1/y(501)^2)*(y(v1-1)-2*y(v1)+y(v1+1));
end
dydt(100)=p.k1_pos*p.clpxp*y(100+200)/(y(300)+p.km1)-p.k1_neg*y(1)+(p.D_complex1/y(501)^2)*(y(99)-y(100));

%%%%%%%%%%%%%%%CpdR_f 

dydt(101)=p.ks_cpdr-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
      +p.k2_pos*y(101+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(101)...
  +p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)*y(101+300)...
      +(p.D_cpdr/y(501)^2)*(y(102)-y(101));
for v2=102:199
  dydt(v2)=p.ks_cpdr-p.kd_cpdr*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
      +p.k2_pos*y(v2+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(v2)...
      +p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2)*y(v2+300)...
      +(p.D_cpdr/y(501)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
end
dydt(200)=p.ks_cpdr-p.kd_cpdr*y(200-100)*y(v2)/(y(200-100)+p.J2)...
      +p.k2_pos*y(200+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(200)...
      +p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)*y(200+300)...
      +(p.D_cpdr/y(501)^2)*(y(199)-y(200));
%%%%%%%%%%%%%%%%CpdR_b
for v3=201:300
  dydt(v3)=-p.kd_cpdr*y(v3-200)*y(v3)/(y(v3-200)+p.J2)...
      +p.kcpdr_f_b*y(v3+200)*y(v3-100)-p.kcpdr_b_f*y(v3)...
  +p.k2_pos*y(v3+100)*p.divkp_oldpole/(p.divkp_oldpole+p.J4)-p.k2_neg*y(v3)...
  +p.k1_neg*y(v3-200)-p.k1_pos*p.clpxp*y(v3)/(y(v3)+p.km1);
end

%%%%%%%%CpdR~P


dydt(301)=p.k2_neg*(y(101)+y(201))-p.k2_pos*y(301)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(302)-y(301))/(y(501)^2);

for v4=302:399
  dydt(v4)=p.k2_neg*(y(v4-200)+y(v4-100))-p.k2_pos*y(v4)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(v4+1)-2*y(v4)+y(v4-1))/(y(501)^2);
end
  dydt(400)=p.k2_neg*(y(400-200)+y(400-100))-p.k2_pos*y(400)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(399)-y(400))/(y(501)^2);

%%CpdR sticky
dydt(401:500)=0;
%%%%%%%%%%%add model 2
%%before uncomment, should change total amount to 701
% %complex2. 
%     dydt(501)=p.k3_pos*y(501-500)*y(501+100)-p.k3_neg*y(501)+p.D_complex2*(y(501+1)-y(501))/(y(701)^2);
% for v5=502:599
%     dydt(v5)=p.k3_pos*y(v5-500)*y(v5+100)-p.k3_neg*y(v5)+p.D_complex2*(y(v5+1)-2*y(v5)+y(v5-1))/(y(701)^2);
% end
%     dydt(600)=p.k3_pos*y(600-500)*y(600+100)-p.k3_neg*y(600)+p.D_complex2*(y(600-1)-y(600))/(y(701)^2);
% %RcdA
%     dydt(601)=p.ks_rcda*y(601)^2/(y(601)^2+p.J5^2)-p.kd_rcda*y(601)+p.D_rcda*(y(601+1)-y(601))/(y(701)^2);
% for v6=602:699
%     dydt(v6)=p.ks_rcda*y(v6)^2/(y(v6)^2+p.J5^2)-p.kd_rcda*y(v6)+p.D_rcda*(y(v6+1)-2*y(v6)+y(v6-1))/(y(701)^2);
% end
%     dydt(700)=p.ks_rcda*y(700)^2/(y(700)^2+p.J5^2)-p.kd_rcda*y(700)+p.D_complex2*(y(700-1)-y(700))/(y(701)^2);

%%%%%%%%%%%%%%%%%%%%%%
dydt(501)=p.growth*y(501);
%%%%


