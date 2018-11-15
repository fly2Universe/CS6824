%odes_div.m


function dydt=odes_div(t,y)
global p;
% T=120;%period of Caulobacter
% t_d=rem(t,T); %return remainder after division t/T
% %p.cckap=p1*t_d.^3 + p2*t_d.^2 + p3*t_d + p4;
p1 =   1.315e-06 ;
       p2 =  -0.0002346 ;
       p3 =    0.001233;
       p4 =      0.9975;
       p.divkp_free=p1*t^3 + p2*t^2 + p3*t + p4;
   %%%%DivK~P_oldpole
    d1 =  -8.515e-06;
       d2 =   0.0009505;
       d3 =    -0.01231 ;
       d4 =      0.1302 ;
     p.divkp_oldpole= (d1*t^3 + d2*t^2 + d3*t + d4).*(t<=98)+0.*(t>98);
  %RcdA
q1 =  -8.725e-08 ;
       q2 =   2.948e-05 ;
       q3 =   -0.003332  ;
       q4 =      0.1384 ;
       q5 =     -0.8569;
       p.rcda=q1*t^4 + q2*t^3 + q3*t^2 + q4*t+ q5;
dydt=zeros(601,1);

%%Complex1 (ClpXP:CpdR)

dydt(1)=p.k1_pos*p.clpxp*y(1+200)/(y(1+200)+p.km1)+p.k1_pos*p.clpxp*y(1+100)/(y(1+100)+p.km1)...
    -2*p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(2)-y(1));
for v1a=2:49
  dydt(v1a)=p.k1_pos*p.clpxp*y(v1a+200)/(y(v1a+200)+p.km1)+p.k1_pos*p.clpxp*y(v1a+100)/(y(v1a+100)+p.km1)...
      -2*p.k1_neg*y(v1a)+(p.D_complex1/y(601)^2)*(y(v1a-1)-2*y(v1a)+y(v1a+1));
end
dydt(50)=p.k1_pos*p.clpxp*y(50+200)/(y(50+200)+p.km1)+p.k1_pos*p.clpxp*y(50+100)/(y(50+100)+p.km1)...
    -2*p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(49)-y(50));
dydt(51)=p.k1_pos*p.clpxp*y(51+200)/(y(51+200)+p.km1)+p.k1_pos*p.clpxp*y(51+100)/(y(51+100)+p.km1)...
    -2*p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(52)-y(51));
for v1b=52:99
  dydt(v1b)=p.k1_pos*p.clpxp*y(v1b+200)/(y(v1b+200)+p.km1)+p.k1_pos*p.clpxp*y(v1b+200)/(y(v1b+200)+p.km1)...
      -2*p.k1_neg*y(v1b)+(p.D_complex1/y(601)^2)*(y(v1b-1)-2*y(v1b)+y(v1b+1));
end
dydt(100)=p.k1_pos*p.clpxp*y(100+200)/(y(100+200)+p.km1)+p.k1_pos*p.clpxp*y(100+200)/(y(100+200)+p.km1)...
    -2*p.k1_neg*y(1)+(p.D_complex1/y(601)^2)*(y(99)-y(100));

%%CpdR_f 
dydt(101)=p.ks_cpdr-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
      +p.k2_pos*y(101+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(101)...
  +p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)*y(101+300)...
      +(p.D_cpdr/y(601)^2)*(y(102)-y(101))...
      +p.k1_neg*y(101-100)-p.k1_pos*p.clpxp*y(101)/(y(101)+p.km1);
for v2a=102:149
  dydt(v2a)=p.ks_cpdr-p.kd_cpdr*y(v2a-100)*y(v2a)/(y(v2a-100)+p.J2)...
      +p.k2_pos*y(v2a+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(v2a)...
      +p.kcpdr_b_f*y(v2a+100)-p.kcpdr_f_b*y(v2a)*y(v2a+300)...
           +p.k1_neg*y(v2a-100)-p.k1_pos*p.clpxp*y(v2a)/(y(v2a)+p.km1)...
      +(p.D_cpdr/y(601)^2)*(y(v2a-1)-2*y(v2a)+y(v2a+1));
end
dydt(150)=p.ks_cpdr-p.kd_cpdr*y(150-100)*y(101)/(y(150-100)+p.J2)...
      +p.k2_pos*y(150+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(150)...
      +p.kcpdr_b_f*y(150+100)-p.kcpdr_f_b*y(150)*y(150+300)...
           +p.k1_neg*y(150-100)-p.k1_pos*p.clpxp*y(150)/(y(150)+p.km1)...
      +(p.D_cpdr/y(601)^2)*(y(149)-y(150));  
dydt(151)=p.ks_cpdr-p.kd_cpdr*y(151-100)*y(101)/(y(151-100)+p.J2)...
      +p.k2_pos*y(151+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(151)...
      +p.kcpdr_b_f*y(151+100)-p.kcpdr_f_b*y(151)*y(151+300)...
           +p.k1_neg*y(151-100)-p.k1_pos*p.clpxp*y(151)/(y(151)+p.km1)...
      +(p.D_cpdr/y(601)^2)*(y(152)-y(151)); 
for v2b=152:199
  dydt(v2b)=p.ks_cpdr-p.kd_cpdr*y(v2b-100)*y(v2b)/(y(v2b-100)+p.J2)...
      +p.k2_pos*y(v2b+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(v2b)...
      +p.kcpdr_b_f*y(v2b+100)-p.kcpdr_f_b*y(v2b)*y(v2b+300)...
           +p.k1_neg*y(v2b-100)-p.k1_pos*p.clpxp*y(v2b)/(y(v2b)+p.km1)...
      +(p.D_cpdr/y(601)^2)*(y(v2b-1)-2*y(v2b)+y(v2b+1));
end  
dydt(200)=p.ks_cpdr-p.kd_cpdr*y(200-100)*y(v2b)/(y(200-100)+p.J2)...
      +p.k2_pos*y(200+200)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(200)...
      +p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)*y(200+300)...
           +p.k1_neg*y(200-100)-p.k1_pos*p.clpxp*y(200)/(y(200)+p.km1)...
      +(p.D_cpdr/y(601)^2)*(y(199)-y(200));

%%CpdR_b
for v3=201:300
  dydt(v3)=-p.kd_cpdr*y(v3-200)*y(v3)/(y(v3-200)+p.J2)...
      +p.kcpdr_f_b*y(v3+200)*y(v3-100)-p.kcpdr_b_f*y(v3)...
  +p.k2_pos*y(v3+100)*p.divkp_oldpole/(p.divkp_oldpole+p.J4)-p.k2_neg*y(v3)...
  +p.k1_neg*y(v3-200)-p.k1_pos*p.clpxp*y(v3)/(y(v3)+p.km1);
end

%%CpdR~P
dydt(301)=p.k2_neg*(y(101)+y(201))-p.k2_pos*y(301)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(302)-y(301))/(y(601)^2);

for v4a=302:349
  dydt(v4a)=p.k2_neg*(y(v4a-200)+y(v4a-100))-p.k2_pos*y(v4a)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(v4a+1)-2*y(v4a)+y(v4a-1))/(y(601)^2);
end
  dydt(350)=p.k2_neg*(y(350-200)+y(350-100))-p.k2_pos*y(350)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(349)-y(350))/(y(601)^2);
dydt(351)=p.k2_neg*(y(151)+y(251))-p.k2_pos*y(351)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(352)-y(351))/(y(601)^2);
for v4b=352:399
  dydt(v4b)=p.k2_neg*(y(v4b-200)+y(v4b-100))-p.k2_pos*y(v4b)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(v4b+1)-2*y(v4b)+y(v4b-1))/(y(601)^2);
end
  dydt(400)=p.k2_neg*(y(400-200)+y(400-100))-p.k2_pos*y(400)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4))...
+p.D_cpdrp*(y(399)-y(400))/(y(601)^2);


%%CpdR sticky
dydt(491:500)=1;
dydt(401:490)=0;

%complex2. 
    dydt(501)=p.k3_pos*y(501-500)*...y(501+100)
        p.rcda-p.k3_neg*y(501)+p.D_complex2*(y(501+1)-y(501))/(y(601)^2);
for v5a=502:549
    dydt(v5a)=p.k3_pos*y(v5a-500)*...y(v5+100)
        p.rcda-p.k3_neg*y(v5a)+p.D_complex2*(y(v5a+1)-2*y(v5a)+y(v5a-1))/(y(601)^2);
end
   dydt(550)=p.k3_pos*y(550-500)*...y(600+100)
        p.rcda-p.k3_neg*y(550)+p.D_complex2*(y(550-1)-y(550))/(y(601)^2);
    dydt(551)=p.k3_pos*y(551-500)*...y(501+100)
        p.rcda-p.k3_neg*y(551)+p.D_complex2*(y(551+1)-y(551))/(y(601)^2);

for v5b=552:599
    dydt(v5b)=p.k3_pos*y(v5b-500)*...y(v5+100)
        p.rcda-p.k3_neg*y(v5b)+p.D_complex2*(y(v5b+1)-2*y(v5b)+y(v5b-1))/(y(601)^2);
end
    dydt(600)=p.k3_pos*y(600-500)*...y(600+100)
        p.rcda-p.k3_neg*y(600)+p.D_complex2*(y(600-1)-y(600))/(y(601)^2);

%%%%%%%%%%%%%%%%%%%%%%
dydt(601)=p.growth*y(601);