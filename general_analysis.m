function [ Soln ] = general_analysis(curr_velo, ball_dim, i, t_cur)

syms t

Soln = [NaN,NaN];

% BEFORE RUNNING THIS, run path.m to initialize variables in the workspace
% at some point, this will be done automatically

%TODO: add friction force components!

%ball_dim = [0.000714,0.0081,1.8738216*10^-8]; %mass,radius, MoI
g = 9.81;

%Ball initial conditions

vgx   = curr_velo(1);
vgy   = curr_velo(2);
omega = curr_velo(3);

%bring in data from path.m
theta_ge = evalin('base','theta');
rad_of_curv_ge = evalin('base','rad_of_curv');
%tmin_ge = evalin('base','tmin'); %Initial t (vector)

%t_cur = tmin_ge(i);

% FOR LOOP START

% FIND t GIVEN x,y

%Calculate agx & agy (need to calc current rad of curv first)
% *** Calc current theta. (0.0001 added since for some cases if t_cur = 0, returns
% error for div by zero)
curr_theta = eval(subs(theta_ge(i),t,t_cur+0.0001));

% *** Calc radius of curvature
try %Can fail to calc RoC in some cases, so check here
    curr_roc = eval(subs(rad_of_curv_ge(i),t,t_cur)); %finds current radius of curvature
catch ME % if you caaaann (feat. Leo DiCaprio)
    curr_roc = 0;
    %disp('Div. by zero was detected during rad_of_curv calculation and fixed. (set to 0)')
end

% *** Calc agy
agy = ((vgx/curr_roc)^2)*curr_roc;
if isnan(agy) %checks if agy is NaN (happens if curr_roc = 0)
    agy = 0;
end

% *** Calc agx
agx = -(ball_dim(1)*g*ball_dim(2)^2*sin(curr_theta))/(ball_dim(3)+ball_dim(1)*ball_dim(2)^2);
Soln(1) = agx;
Soln(2) = agy;
Soln(3) = curr_theta;
end




