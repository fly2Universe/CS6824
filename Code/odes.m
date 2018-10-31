%odes.m


function dydt=odes(t,y);
global p;

dydt=zeros(301,1);

%%Complex1 (ClpXP:CpdR)
dydt(1)=p.k1_pos*p.clpxp*y(101)-p.k1_neg*y(1);
for v1=2:99
  dydt(v1)=p.k1_pos*p.clpxp*y(v1+100)-p.k1_neg*y(v1);
end
dydt(100)=p.k1_pos*p.clpxp*y(200)-p.k1_neg*y(100);

%%CpdR
for v2=101:200
  dydt(v2)=p.ks_cpdr*y(v2+100)/(y(v2+100)+p.J1)+p.k2_pos*y(v2+100)-p.kd*y(v2-100)*y(v2)/(y(v2-100)+p.J2)-p.k2_neg*y(v2+100)*y(v2)/(y(v2+100)+p.J3);
end

%%CpdR~P
for v3=201:300
  dydt(v3)=p.k2_neg*y(v3)*y(v3-100)/(y(v3)+p.J3)-p,k2_pos*y(v3);
end

%%%%%%%%%%%%%%%%%%%%%%
dydt(301)=p.growth*y(301);