
clf
clear all

path_param = NaN([1,9]);
    path_param(1) = 0.075; %start_hieght
    path_param(2) = 1.25; %ramp_angle 
    path_param(3) = 0.02; %ramp_length 
    path_param(4) = 2.5; %arc_1_dia 
    path_param(5) = 0.085; %landing_height
    path_param(6) = 5; %arc_2_dia
    path_param(7) = 0.025; %flat_legnth 
    path_param(8) = 7; %segments
    path_param(9) = 0.001; %tolerance
    
ball_dim = NaN([1,3]);
    ball_dim(1) = 0.000551; %mass (kg)
    ball_dim(2) = 0.0081; %raduis (m)
    ball_dim(3) = (7/5)*ball_dim(1)*ball_dim(2)^2; %monemt of inertia
    
anal_param = NaN([1,2]);
    anal_param(1) = 0.005; %step size in seconds
    anal_param(2) = 0.00001; %atan error corection
    anal_param(3) = 0.017; %coeifiecnt of friction
    
path_array = path(path_param,ball_dim);
rad_of_curv = rad_of_curv_funct(path_array(1,:),path_array(2,:));
theta = theta_funct(path_array(1,:),path_array(2,:));   

kinn_anal(anal_param,ball_dim,path_array,theta,rad_of_curv)
    
