%Function to test linear_analysis of short linear section ***in reverse!***
x_start = 10 - (75*sin(157/40))/2 + 13207/80;
x_end =   0 - (75*sin(157/40))/2 + 13207/80;

y_start = (75*cos(157/40))/2 + (10*sin(157/40))/(cos(157/40) - 1) + 844;
y_end =   (75*cos(157/40))/2 + (0*sin(157/40))/(cos(157/40) - 1) + 844;
%****************************so many constants



%x = 195;
step = 0.01; %in seconds
ball_dim = [0.0007145489,0.0081,1.87526213316*10^-8];

%find initial conditions (assume velocity = 0)
%x = x_start; 
v = 0; %initial velocity
%solve for y-position

syms t;
equn = x == t - (75*sin(157/40))/2 + 13207/80;
ti = solve(equn,t,'Real',1);
y = (75*cos(157/40))/2 + (ti*sin(157/40))/(cos(157/40) - 1) + 844;
%Calculate deriv & theta (constant in this case)
deriv = sin(157/40)/(cos(157/40) - 1);
%theta = atan(deriv);
theta = -8.2463*2*pi/360

x = 0;
x_end = 700;
x(1)=0
iter = 2;

while x < x_end
    %Find accel of ball
    A = linear_analysis(theta, ball_dim); %vector containing accel, ang_accel data
    %calculate new velocity & position
    v_new = v + A(1)*step;
    v = v_new;
    
    %over very short distances & periods of time, we can assume
    %acceleration is constant (irl, it will vary)
    %for this case, accel is constant as it is a linear slope
    x_new = x(iter-1) + v_new*step + 0.5*A(1)*step^2;
    x(iter) = x_new;
    y(iter) = 100 + x*(10/69);
    iter = iter + 1
end
x 
v
iter*step
for iter2 = 1:length(x)
    hold on
    plot(x(iter2),y(iter2));
end

%r_curv = 