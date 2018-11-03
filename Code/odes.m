%odes.m


function dydt=odes(t,y);
global p;

dydt=zeros(301,1);

%%Complex1 (ClpXP:CpdR)
dydt(1)=p.k1_pos*p.clpxp*y(1+200)-p.k1_neg*y(1)+(p.D_complex1/y(501)^2)*(y(2)-y(1));
for v1=2:99
  dydt(v1)=p.k1_pos*p.clpxp*y(v1+200)-p.k1_neg*y(v1)+(p.D_complex1/y(501)^2)*(y(v1-1)-2*y(v1)+y(v1+1));
end
dydt(100)=p.k1_pos*p.clpxp*y(100+200)-p.k1_neg*y(1)+(p.D_complex1/y(501)^2)*(y(99)-y(100));

%%CpdR_f 
dydt(101)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(101-100)*y(101)/(y(101-100)+p.J2)...
      +p.k2_pos*y(101+200)-p.k2_neg*p.cckap*y(101)/(p.cckap+p.J3)+p.kcpdr_b_f*y(101+100)-p.kcpdr_f_b*y(101)...
      +(p.D_cpdrf/y(501)^2)*(y(102)-y(101));
for v2=102:199
  dydt(v2)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
      +p.k2_pos*y(v2+200)-p.k2_neg*p.cckap*y(v2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2)...
      +(p.D_cpdrf/y(501)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
end
dydt(200)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd_cpdr*y(200-100)*y(v2)/(y(200-100)+p.J2)...
      +p.k2_pos*y(200+200)-p.k2_neg*p.cckap*y(200)/(p.cckap+p.J3)+p.kcpdr_b_f*y(200+100)-p.kcpdr_f_b*y(200)...
      +(p.D_cpdrf/y(501)^2)*(y(199)-y(200));
%%CpdR_b
for v3=201:300
  dydt(v3)=-p.kd_cpdr*y(v3-200)*y(v3)/(y(v3-200)+p.J2)+p.kcpdr_f_b*y(v3-100)-p.kcpdr_b_f*y(v3)...
      +p.k2_pos*y(v3+200)-p.k2_neg*p.cckap*y(v3)/(p.cckap+p.J3)+p.k1_neg*y(v3-200)-p.k1_pos*y(v3)*p.clpxp;
end


%%CpdR~P_f
dydt(301)=p.k2_neg*p.cckap*y(101)/(p.cckap+p.J3)-p.k2_pos*y(301)+p.kcpdrp_b_f*y(401)-p.kcpdrp_f_b*y(301)...
+p.D_cpdrp*(y(302)-y(301))/(y(501)^2);

for v4=302:399
  dydt(v4)=p.k2_neg*p.cckap*y(v4-200)/(p.cckap+p.J3)-p.k2_pos*y(v4)+p.kcpdrp_b_f*y(v4+100)-p.kcpdrp_f_b*y(v4)...
+p.D_cpdrp*(y(v4+1)-2*y(v4)+y(v4-1))/(y(501)^2);
end
dydt(400)=p.k2_neg*p.cckap*y(200)/(p.cckap+p.J3)-p.k2_pos*y(400)+p.kcpdrp_b_f*y(500)-p.kcpdrp_f_b*y(400)...
+p.D_cpdrp*(y(399)-y(400))/(y(501)^2);


%%CpdR~P_b
for v5=401:500
  dydt(v5)=p.k2_neg*p.cckap*y(v5-200)/(p.cckap+p.J3)-p.k2_pos*y(v5)+p.kcpdrp_f_b*y(v5-100)-p.kcpdrp_b_f*y(v5);
end


%%%%%%%%%%%%%%%%%%%%%%
dydt(501)=p.growth*y(501);