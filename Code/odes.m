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
  dydt(v2)=p.ks_cpdr*y(v2+100)/(y(v2+100)+p.J1)+p.k2_pos*y(v2+100)-p.kd*y(v2-100)*y(v2)/(y(v2-100)+p.J2)-p.k2_neg*y(v2+100)*y(v2)/(y(v2+100)+p.J3);
end

%%CpdR_b
for v3=201:300
  dydt(v3)=


%%CpdR~P_f
for v4=301:400
  dydt(v4)=p.k2_neg*p.cckap*y(v4-200)/(p.cckap+p.J3)-p.k2_pos*y(v4)+p.kcpdrp_b_f*y(v4+100)-p.kcpdrp_f_b*y(v4);
end


%%CpdR~P_b
for v5=401:500
  dydt(v5)=p.k2_neg*p.cckap*y(v5-200)/(p.cckap+p.J3)-p.k2_pos*y(v5)+p.kcpdrp_f_b*y(v5-100)-p.kcpdrp_b_f*y(v5);
end

%%%%%%%%%%%%%%%%%%%%%%
dydt(301)=p.growth*y(301);