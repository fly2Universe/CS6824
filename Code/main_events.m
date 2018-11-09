%main_events.m

function sw_to_early_pd = main_events(initial_cond)
param(1); %cell can grow

%% INITIAL CONDITIONS
y0 = initial_cond.init;
y0 = y0.';
%% INTEGRATION PARAMETERS
t0 = 0;
tf = 120;

options = odeset('Events',@event,'RelTol',1e-4,'AbsTol',1e-6);
tout=t0;
yout=y0.';
teout = [];
yeout = [];
ieout = [];

while t0<tf
    
    [t,y,te,ye,ie]=ode15s(@odes,[t0 tf],y0,options);
    
    nt=length(t);
    
    tout=[tout;t(2:nt)];
    yout=[yout;y(2:nt,:)];
    teout = [teout;te];
    yeout = [yeout;ye];
    ieout = [ieout;ie];
    
    y0 = y(nt,:);

    if isscalar(ie) == 0
        ie = 0;
    end
    
   
    %enforced localization Sticky CpdR
    if (ie==1)||(ie==2)
        y0(501:511)=1;
        y0(512:600)=0;
    elseif ie==3
        y0(501:600)=0;
    end    
    
    t0=t(nt);
    if t0>=tf
        break;
    end
    
% Define Grid M
for n=51:100
    M(:,n)=yout(:,4901)*(n/100)-0.5*(yout(:,4901)*(.5) + yout(:,4901)*.51); 
    % Each element from column 51 to 100 is assigned the 'n'th fraction of 
    %total cell length at a given time step followed by subtracting the mid 
    %point value so that the centre of the grid takes the value 0
end

M(:,1:50)=fliplr(M(:,51:100));      % The 1 st 50 colums is the flipped version of the next 50 eg: 3 2 1 1 2 3
M(:,1:50)=-M(:,1:50);               % now its -3 -2 -1 1 2 3
M=100*M;



figure(1)
ax1 = subplot(2,2,1);
pcolor(tout, M, plec_kinase)
shading interp
colorbar
title('PleC kinase')
xlim([0 120])
label_str = strcat('cell size (',char(956),'m)');
ylabel(label_str) 

ax2 = subplot(2,2,2);
pcolor(tout, M, divk_p)
shading interp
colorbar
xlim([0 120])
title('DivK~P')

ax3 = subplot(2,2,3);
pcolor(tout, M, divl_free)
shading interp
colorbar
xlim([0 120])
xlabel('time (min)')
title('active DivL')
ylabel(label_str) 

ax4 = subplot(2,2,4);
pcolor(tout, M, ctra_p)
shading interp
colorbar
xlim([0 120])
xlabel('time (min)')
title('CtrA~P')

h = suptitle(initial_cond.strain);
set(h,'interpreter','none')


    sw_to_early_pd.endpoints = yout(end,:);

end