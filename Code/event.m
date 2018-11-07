function [value,isterminal,direction] = event(t,y)

sw_end = 30;%50;
st_start = 50;%%70;
pd = 90;%100;
dummy = 200;

value=[sign(t-sw_end); sign(t-st_start); sign(t-pd); sign(t-dummy)];
isterminal=[1;1;1;1];
direction=[+1;+1;+1;+1];
end