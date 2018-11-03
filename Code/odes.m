%odes.m


function dydt=odes(t,y);
global p;

dydt=zeros(301,1);

%%Complex1 (ClpXP:CpdR)
for v1=1:100
  dydt(v1)=p.k1_pos*p.clpxp*y(v1+200)-p.k1_neg*y(v1);
end

%%CpdR_f 
for v2=101:200
  dydt(v2)=p.ks_cpdr*p.cckap/(p.cckap+p.J1)-p.kd*y(v2-100)*y(v2)/(y(v2-100)+p.J2)...
      +p.k2_pos*y(v2+200)-p.k2_neg*p.cckap*y(v2)/(p.cckap+p.J3)+p.kcpdr_b_f*y(v2+100)-p.kcpdr_f_b*y(v2);
end

%%CpdR_b
for v3=201:300
  dydt(v3)=-p.kd*y(v3-200)*y(v3)/(y(v3-200)+p.J2)+p.kcpdr_f_b*y(v3-100)-p.kcpdr_b_f*y(v3)...
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