clc
clear
clear all

%TODO Change to Energy Equn method

%Doesn't work at the moment, will change to energy equn method and try
%again

m = 1; %In kg
g = 9.81;
theta = 15; %deg, can change to rad later
radius = 0.01; %In meters
Ig = 4*10^-7

syms agx agy %Setup of agx and agy unknown variables for symbolic equation

eqns = [-m*g*radius*sind(theta) == -(radius*cosd(theta))*m*agx + (radius*sind(theta))*m*agy + Ig*(sqrt(agx^2+agy^2))/radius, cosd(theta)*agx == sind(theta)*agy]
S = solve(eqns, [agx,agy])
S.agx
S.agy