
clf
clear all

parameters = NaN([1,9]);
    parameters(1) = 0.075; %start_hieght
    parameters(2) = 1.25; %ramp_angle 
    parameters(3) = 0.02; %ramp_length 
    parameters(4) = 2.5; %arc_1_dia 
    parameters(5) = 0.085; %landing_height
    parameters(6) = 5; %arc_2_dia
    parameters(7) = 0.025; %flat_legnth 
    parameters(8) = 7; %segments
    parameters(9) = 0.001; %tolerance
    
ball_dim = NaN([1,3]);
    ball_dim(1) = 0.000551; %mass (kg)
    ball_dim(2) = 0.0081; %raduis (m)
    ball_dim(3) = (7/5)*ball_dim(1)*ball_dim(2)^2; %monemt of inertia
    
path_array = path(parameters,ball_dim);
    
    
