function [ Soln ] = linear_analysis(theta, ball_dim)
%Code that finds acceleration and angular accel of ball given theta and
%ball dimensions

%Input ball dimensions (for ease of reading)
m = ball_dim(1); %mass
radius = ball_dim(2);
Ig = ball_dim(3);

g = 9.81;
%angle theta is given as input


% Below was used for testing function VVV
%m = 1; %In kg
%g = 9.81;
%theta = pi/2; %in radians (positive = angled up like this:  /  )
%radius = 0.01; %In meters
%Ig = 4*10^-7; %guesstimated


syms apx %Setup of apx unknown variable for symbolic equation

equn = -m*g*sin(theta)*radius == +radius*m*apx + (Ig+m*(radius^2))*(apx/radius);

Soln(1) = solve(equn, apx);

Soln(2) = -Soln(1)/radius; %gives ang_accel direction according to right-hand rule

end

%don't use this equation
%eqns = [-m*g*radius*sind(theta) == -(radius*cosd(theta))*m*agx + (radius*sind(theta))*m*agy + Ig*(sqrt(agx^2+agy^2))/radius, cosd(theta)*agx == sind(theta)*agy]
