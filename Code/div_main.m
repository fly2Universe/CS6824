unction early_to_late_pd = div_main(sw_to_early_pd) 

param(1, sw_to_early_pd.type)

%% INITIAL CONDITIONS

y0 = sw_to_early_pd.endpoints;
y0 = y0.';
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 30;

tout=t0;
y=y0;
     
[t,y]=ode15s(@odes_div,[t0 tf],y0);

early_to_late_pd.endpoints = y(end,:);

%% Plots
for n=51:100
    M(:,n)=y(:,4901)*(n/100)-0.5*(y(:,4901)*(.5) + y(:,4901)*.51); % Each element from column 51 to 100 is assigned the 'n'th fraction of total cell length at a given time step followed by subtracting the mid point value so that the centre of the grid takes the value 0
end

M(:,1:50)=fliplr(M(:,51:100)); % The 1 st 50 colums is the flipped version of the next 50 eg: 3 2 1 1 2 3
M(:,1:50)=-M(:,1:50); %now its -3 -2 -1 1 2 3
%  The Y-axis variables are recorded in matrix G

M=100*M;

plec_kinase(:,1:100)=y(:,601:700)+y(:,701:800)+y(:,801:900)+y(:,901:1000)+y(:,1001:1100)+y(:,1101:1200)+y(:,1201:1300)+y(:,1301:1400)+y(:,1401:1500)+y(:,1501:1600)+y(:,1601:1700)+y(:,1701:1800)+y(:,1801:1900)+y(:,1901:2000)+y(:,2101:2200);
divk_p(:,1:100)=y(:,1:100)+y(:,501:600)+2*y(:,601:700)+2*y(:,701:800)+y(:,801:900)+y(:,901:1000)+y(:,1001:1100)+2*y(:,1101:1200)+y(:,1201:1300)+y(:,1501:1600)+y(:,1901:2000)+y(:,2601:2700) + y(:,2901:3000);
ccka_kinase(:,1:100)=y(:,3501:3600);%+y(:,3601:3700)+y(:,3701:3800)+y(:,4101:4200)+y(:,4201:4300);
ctra_p(:,1:100) = y(:,3201:3300);

plec_kinase = fliplr(plec_kinase);
ccka_kinase = fliplr(ccka_kinase);
divk_p = fliplr(divk_p);
ctra_p = fliplr(ctra_p);

plec_kinase = plec_kinase.';
ccka_kinase = ccka_kinase.';
divk_p = divk_p.';
ctra_p = ctra_p.';
M = M.';

figure(2)
ax1 = subplot(2,2,1);
pcolor(t, M, plec_kinase)
shading interp
colorbar
title('PleC kinase')
xlim([0 30])
label_str = strcat('cell size (',char(956),'m)');
ylabel(label_str) 

ax2 = subplot(2,2,2);
pcolor(t, M, divk_p)
shading interp
colorbar
xlim([0 30])
title('DivK~P')

ax3 = subplot(2,2,3);
pcolor(t, M, ccka_kinase)
shading interp
colorbar
xlim([0 30])
xlabel('time (min)')
title('CckA kinase')
ylabel(label_str) 

ax4 = subplot(2,2,4);
pcolor(t, M, ctra_p)
shading interp
colorbar
xlim([0 30])
xlabel('time (min)')
title('CtrA~P')

h = suptitle(sw_to_early_pd.strain);
set(h,'interpreter','none')
end