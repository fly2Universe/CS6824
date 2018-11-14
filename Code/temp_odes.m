%temperal_odes.m

function dydt=temp_odes(t,y)
global p;

 T=120;%period of Caulobacter before division
 t_d=rem(t,T); %return remainder after division t/T
%%%%DivK~P_free
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
 %%%%%%Divk 120min
 q1 =  -1.843e-06;
       q2 =   0.0003078 ;
       q3 =    -0.01197 ;
       q4 =      0.9516 ;
          p.divkp_120=q1*t_d.^3 + q2*t_d.^2 + q3*t_d + q4;
%     a1 =     -0.1011 ;
%        c1 =     0.07546;
%        d1=      0.9291 ;
% p.divkp_120=a1*sin(3*pi*t_d./180+c1)+d1;
%RcdA
 q1 =  -8.725e-08 ;
       q2 =   2.948e-05 ;
       q3 =   -0.003332  ;
       q4 =      0.1384 ;
       q5 =     -0.8569;
       p.rcda=q1*t_d.^4 + q2*t_d.^3 + q3*t_d.^2 + q4*t_d+ q5;
%%Complex1 (ClpXP:CpdR)
  dydt(1)=p.k1_pos*p.clpxp*y(3)/(y(3)+p.km1)-p.k1_neg*y(1);

%%CpdR_f 
dydt(2)=p.ks_cpdr-p.kd_cpdr*y(1)*y(2)/(y(1)+p.J2)+p.k2_pos*y(4)*p.divkp_free/(p.divkp_free+p.J3)-p.k2_neg*y(2)...
    +p.kcpdr_b_f*y(3)-p.kcpdr_f_b*y(2);

%%CpdR_b
  dydt(3)=-p.kd_cpdr*y(1)*y(3)/(y(1)+p.J2)+p.kcpdr_f_b*y(2)-p.kcpdr_b_f*y(3)...
      +p.k2_pos*y(4)*p.divkp_oldpole/(p.divkp_oldpole+p.J4)-p.k2_neg*y(3)...
      +p.k1_neg*y(1)-p.k1_pos*p.clpxp*y(3)/(y(3)+p.km1);
%%%%%%use divk 120
% dydt(2)=p.ks_cpdr-p.kd_cpdr*y(1)*y(2)/(y(1)+p.J2)+p.k2_pos*y(3)*p.divkp_120/(p.divkp_120+p.J3)-p.k2_neg*y(2);
%%CpdR~P
  dydt(4)=p.k2_neg*(y(2)+y(3))-p.k2_pos*y(4)*(p.divkp_free/(p.divkp_free+p.J3)+p.divkp_oldpole/(p.divkp_oldpole+p.J4));
%%%%use divk 120
 % dydt(3)=p.k2_neg*y(2)-p.k2_pos*y(3)*p.divkp_120/(p.divkp_120+p.J3);

%%%%%%%%%%%add model 2
%complex2
%dydt(5)=p.k3_pos*y(1)*y(6)-p.k3_neg*y(5);

dydt(5)=p.k3_pos*y(1)*p.rcda-p.k3_neg*y(5);
% %RcdA
% dydt(6)=p.ks_rcda*y(6)^2/(y(6)^2+p.J5)-p.kd_rcda*y(6)*y(1)/(y(1)+p.J6);

%%%%%%%%%%%%%%%%%%%%%%
%%%%
dydt=dydt';

