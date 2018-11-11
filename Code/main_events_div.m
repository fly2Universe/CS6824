%main_events_div.m

function main_events_div(sw_to_early_pd)
param(1); %cell can grow

%% INITIAL CONDITIONS
y0 = sw_to_early_pd;
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


%simulate the system in a segmented fashion, which includes a series of
%event
while t0<tf
    %solve system, integration is stopped when a specific event happens
    [t,y,te,ye,ie]=ode15s(@odes_div,[t0 tf],y0,options);
    %number of points calculated
    nt=length(t);
    %the time points and correspoding values calculated before the event
    tout=[tout;t(2:nt)];
    yout=[yout;y(2:nt,:)];
    %record the time when an event happens, the correspoding y value and the event index 
    teout = [teout;te];
    yeout = [yeout;ye];
    ieout = [ieout;ie];
    
    %reset the inition conditon
    y0 = y(nt,:);
    
    %if something wrong
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

end


% Define Grid M for plot
% The plot is basically a heat map of the concentration of the species along the main axis of the cell
for n=51:100
    M(:,n)=yout(:,601)*(n/100)-0.5*(yout(:,601)*(.5) + yout(:,601)*.51); 
    % Each element from column 51 to 100 is assigned the 'n'th fraction of 
    %total cell length at a given time step followed by subtracting the mid 
    %point value so that the centre of the grid takes the value 0
end

M(:,1:50)=fliplr(M(:,51:100));      % The 1 st 50 colums is the flipped version of the next 50 eg: 3 2 1 1 2 3
M(:,1:50)=-M(:,1:50);               % now its -3 -2 -1 1 2 3
M=100*M;

%calculation of concentration
%complex1
complex1(:,1:100)=yout(:,1:100);
%free+bound CpdR
CpdR(:,1:100)=yout(:,101:200)+yout(:,201:300);
%free+bound CpdR~p
CpdR_p(:,1:100)=yout(:,301:400)+yout(:,401:500);

%flip array ?
complex1 = fliplr(complex1);
CpdR = fliplr(CpdR);
CpdR_p = fliplr(CpdR_p);

%transpose ?
complex1 = complex1.';
CpdR = CpdR.';
CpdR_p = CpdR_p.';
M = M.';


figure(1)
ax1 = subplot(2,2,1);
pcolor(tout, M, complex1)
shading interp
colorbar
title('ClpXP:CpdR')
xlim([0 120])
xlabel('time (min)')
label_str = strcat('cell size (',char(956),'m)');
ylabel(label_str) 

ax2 = subplot(2,2,2);
pcolor(tout, M, CpdR)
shading interp
colorbar
xlim([0 120])
xlabel('time (min)')
label_str = strcat('cell size (',char(956),'m)');
ylabel(label_str) 
title('CpdR')

ax3 = subplot(2,2,3);
pcolor(tout, M, CpdR_p)
shading interp
colorbar
xlim([0 120])
xlabel('time (min)')
title('CpdR~p')
label_str = strcat('cell size (',char(956),'m)');
ylabel(label_str) 







end