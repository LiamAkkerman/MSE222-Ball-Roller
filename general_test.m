% general_test.m

% get results for slowest, 6 for each parameter, fastest, regular, mean and
% find time for each one


%TODO: Add friction components (can be done in general_analysis?)
% I think that when the theta is evaluated, it kinda goes a bit crazy when
% it's approaching vertical angles like 90deg or -90deg which is causing
% the weird spikes on the plot and on the accel graphs?


% ############
% This code determines the time it takes for a marble of defined dimensions
% to travel along all 7 predefined segments, as well as determining the
% velocity, accel, and angular velo/accel at each iteration point.

% This is essentially done by determining the acceleration of the ball at a
% given point, then by assuming that the acceleration would be "constant"
% over a very short time interval, finding the resulting new velocity and
% ball position.

% For the actual calculations, the acceleration is found by finding the
% current t_curr value for that point on the line. Basically, the lines
% representing the curves are parametric equations, so t_curr is the value
% that represents the current X & Y coordinates of the point of contact of
% the marble. Knowing t_curr also gives us the normal and tangential
% acceleration of the ball as well as some other useful stuff too. To find
% the location of the center of mass of the marble, we can take the
% position of the point of contact and add the offset from the radius of
% the marble if we know theta (can be found using t_curr). By using
% acceleration at that point, we can find the new X Y coords of the center
% of mass of the marble. We then need to find the respective t_curr value
% for the new point of contact of the marble with the line. We can find the
% approximate location of that point by subtracting the offset of the radius of
% the ball by knowing the previous angle plus a guesstimated value
% 'theta_change' to improve accuracy (not used in segment 4 atm). Now that
% the approximate location of the point of contact is found, we can find
% the t_curr value by using the equation (X_contact = sx_ge(#)),
% essentially meaning we can use the x-component of the parametric equation
% and the x value of the contact point to find t_curr (as long as the 
% segment passes the vertical line test, which they all do). Now that we have 
% t_curr, the loop checks if it's reached the end of the segment and the 
% process continues. :)

% Note that in the calculations, lowercase (x,y) represent tangential and 
% normal components and uppercase (X,Y) represent global coordinates. Ex:
% vgx_new is the tangential velocity of the ball. Sometimes the coordinate
% systems are weird for some segments so be careful lol
% ############

clc
clear all

path()


%Setup initial conditions
curr_velo = [0,0,0,0]; %vgx, vgy, omega, alpha

%Normal Parameters
ball_dim = [0.000551,0.0081,NaN,0.017]; %mass (kg), radius (m), MoI (kg/m^2), coeff of friction
%ball_dim(4) = ball_dim(4)*1.1;
%(Parameters * 90%)
%ball_dim = [0.000551,0.0081,NaN,0.017]*0.9; %mass (kg), radius (m), MoI (kg/m^2), coeff of friction

%Parameters *110%
ball_dim = [0.000551,0.0081,NaN,0.017]*0.9; %mass (kg), radius (m), MoI (kg/m^2), coeff of friction



ball_dim(3) = (2/5)*ball_dim(1)*(ball_dim(2))^2;
timePerIter = 0.002;  %step size (in seconds)
iter = 1;
syms tvar;
g = 9.81;

%Bring in info from path.m
tmin_ge = evalin('base','tmin'); %Initial t (vector) for each segment
tmax_ge = evalin('base','tmax'); %final t (vector) for each segment
t_curr = tmin_ge(1); %Starts at segment 1, so find initial value for 1st seg

sx_ge = evalin('base','sx'); %x-parametric equation
sy_ge = evalin('base','sy'); %y-parametric equation
theta_ge = evalin('base','theta');

X_curr = eval(subs(sx_ge,t,t_curr));
X_curr = X_curr(1) + ball_dim(2)*cos(eval(subs(theta_ge(1),t,t_curr+0.001))+pi/2); %Finds initial x,y coords of the center of mass
Y_curr = eval(subs(sy_ge,t,t_curr));
Y_curr = Y_curr(1) + ball_dim(2)*sin(eval(subs(theta_ge(1),t,t_curr+0.001))+pi/2);
%this is found by finding the position of contact on the curve, then adding
%the radius of the ball to find the CoM (radius is in m, so convert to mm)
theta_prev = eval(subs(theta_ge(1),t,t_curr+0.001)); %used later for accuracy improvement for finding X_contact




% NOTE: lowercase x,y coordinates represent the tangential and normal
% components of the ball for velocity, acceleration, etc analysis
% UPPERCASE X,Y coordinates represent the global X Y coords (X = left to
% right, Y = down to up) for position analysis!

% ***** First segment
while t_curr < tmax_ge(1)
    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 1,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (x)
    %vgy_new = curr_velo(2) + a_new(2)*timePerIter; %normal (y) NOTE: this
    %part is incorrect! there is no instantaneous normal velocity of
    %the ball, but there IS normal acceleration!
    vgy_new = 0; % No normal component of velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    %***Find X_new, Y_new
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new; %redundant, can be removed?
    %This is done by converting vgx & vgy tangential and normal velocities to the horizontal &
    %vertical coord system to find X_new and Y_new.
    
    %***Find change in theta between previous iterations (used to improve 
    %accuracy of finding t_curr)
    theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_prev = a_new(3);
    
    %***Find new 't_curr' by calculating the x-coord of the point of
    %contact
    X_contact = X_curr - ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %Some weird geometry stuff going on here. Basically, this accounts for
    %how the change of the position of the center of mass isn't exactly the
    %same as the change of the position of the point of contact. Pos. of
    %the contact point is used for finding the next t_curr value, so it's
    %fairly important to account for the difference.
    
    t_curr = eval(solve(-X_contact + (3*tvar)/80 - (3*sin(tvar))/80 + 173/10000)) %uses sx_ge(1) equation
    
    %*** 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(1),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
   
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    %Old way of recording values
%     results(1,iter) = X_curr; % x position using global coords
%     results(2,iter) = Y_curr; % y position using global coords
%     results(3,iter) = sqrt(curr_velo(1)^2+curr_velo(2)^2); % magnitude of velo
%     results(4,iter) = sqrt(a_new(1)^2 + a_new(2)^2); % magnitude of accel
%     results(5,iter) = curr_velo(3); % angular velocity
%     results(6,iter) = curr_velo(4); % angular acceleration
%     results(7,iter) = a_new(3); %theta
%     results(8,iter) = a_new(1);
%     results(9,iter) = a_new(2);
   
    
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter


% ***** Second Segment
t_curr = tmin_ge(2); %this isn't always true for the start of each segment!

while t_curr < tmax_ge(2)

    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 2,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (lowercase x)
    vgy_new = 0; % no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    
    %***Find x_new, y_new (note: a_new(3) is theta)
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new;
    %This is done by converting vgx & vgy tangential and normal velocities
    %(they should probably called something else lol) to the horizontal &
    %vertical coord system to find X_new and Y_new. X & Y are in mm, so
    %convert to meters first.
    
    %***Find change in theta between previous iterations (used to improve 
    %accuracy of finding t_curr)
    theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_prev = a_new(3);
    
    X_contact = X_curr - ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %For explination, see segment 1 description of this part
    
    %*************************************again, what are these numbers?
    t_curr = X_contact - 34398703330718350307/180143985094819840000 %uses sx_ge(2) equation
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(2),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(1,iter) = X_curr; % x position using global coords
    results(2,iter) = Y_curr; % y position using global coords
    results(3,iter) = sqrt(curr_velo(1)^2+curr_velo(2)^2); % magnitude of velo
    results(4,iter) = sqrt(a_new(1)^2 + a_new(2)^2); % magnitude of accel
    results(5,iter) = curr_velo(3); % angular velocity
    results(6,iter) = curr_velo(4); % angular acceleration
    results(7,iter) = a_new(3); %theta
    
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter


% ***** Segment 3 (from here on, minimal comments to reduce visual clutter)
t_curr = tmin_ge(3)

while t_curr < tmax_ge(3)

    % *** Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 3,t_curr); %finds agx,agy,theta
    
    % *** Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (x)
    vgy_new = 0; %no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    % *** Find X_new, Y_new
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new; 

    % *** Theta Change Calculation
    % In this section, ignore theta_change correction since it makes the
    % whole thing a bit wonky for now ? fix later?
    %theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_change = 0;
    theta_prev = a_new(3);
    
    %*** Find new 't_curr'
    % Completed by calculating the x-coord of the point of contact
    X_contact = X_curr + ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %Weird geometry going on here, see segment 1 or 2 for info. Note, this
    %is changed compared to seg 1 & 2 versions due to track inversion,
    
    % *** Calculate t_curr
    t_curr = acos(321321996768819125/30399297484750848 - (4000*X_contact)/81)
    %Above equation was found by solving the sx_ge(3) == X_contact equation for t 
    
    % *** Snap Y_curr to line to account for slight drift 
    Y_curr = eval(subs(sy_ge(3),t,t_curr)) + ball_dim(2)*sin(a_new(3)-pi/2);
   
    % *** Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1; 
end
time_taken = (iter-1)*timePerIter


% ***** Fourth Segment
t_curr = tmin_ge(4); %this isn't always true for the start of each segment!

curr_velo(1) = -curr_velo(1); % To account for ball travelling right to left
theta_prev = 1.5707; % To account for theta_change being unused in previous segment

while t_curr < tmax_ge(4)

    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 4,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (lowercase x)
    vgy_new = 0; % no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    
    %***Find x_new, y_new (note: a_new(3) is theta)
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new;
    
    %***Find change in theta between previous iterations
    theta_change = a_new(3) - theta_prev; 
    theta_change = 0;
    theta_prev = a_new(3);
    
    X_contact = X_curr - ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %For explination, see segment 1 description of this part
    
    t_curr = solve ((17*sin(t))/400 - (17*t)/400 + 337652442483427206561/1441151880758558720000-X_contact,t)
    
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(4),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter


% ***** Segment 5

t_curr = tmax_ge(5)

%NOTE!!: For this segment, t goes from tmax to tmin for some reason!

while t_curr > tmin_ge(5)

    % *** Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 5,t_curr); %finds agx,agy,theta
    
    % *** Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (x)
    vgy_new = 0; %no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    % *** Find X_new, Y_new
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new; 

    % *** Theta Change Calculation
    % In this section, ignore theta_change correction since it makes the
    % whole thing a bit wonky for now ? fix later?
    %theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_change = 0;
    theta_prev = a_new(3);
    
    %*** Find new 't_curr'
    % Completed by calculating the x-coord of the point of contact
    X_contact = X_curr + ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %Weird geometry going on here, see segment 1 or 2 for info. Note, this
    %is changed compared to seg 1 & 2 versions due to track inversion,
    
    % *** Calculate t_curr
    t_curr = acos(908926702018138375/364791569817010176 - (2000*X_contact)/81);
    %Above equation was found by solving the sx_ge(5) == X_contact equation for t 
    
    % *** Snap Y_curr to line to account for slight drift 
    Y_curr = eval(subs(sy_ge(5),t,t_curr)) + ball_dim(2)*sin(a_new(3)-pi/2);
   
    % *** Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1; 
end



% ***** Segment 6
t_curr = tmin_ge(6); %this isn't always true for the start of each segment!

curr_velo(1) = -curr_velo(1); % To account for ball travelling left to right
theta_prev = -1.5704;

while t_curr < tmax_ge(6)

    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 6,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (lowercase x)
    vgy_new = 0; % no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    
    %***Find x_new, y_new (note: a_new(3) is theta)
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new;
    
    %***Find change in theta between previous iterations 
    theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_change = 0;
    theta_prev = a_new(3);
    
    X_contact = X_curr - ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %For explination, see segment 1 description of this part
    
 
    %Find t_curr using sx_ge(6) equation (simplified) and X_contact
    t_curr = solve((11350419166481883404180349057063*t)/44328935437225859658807934713856 - (11350419166481883404180349057063*sin(t))/44328935437225859658807934713856 + 5441351322011282287/90071992547409920000-X_contact)
    
    
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(6),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter

% ***** Segment 7
t_curr = tmin_ge(7); %this isn't always true for the start of each segment!

while X_curr < 0.889 % While ball position isn't at the finish position

    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 7,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*timePerIter; %tangential (lowercase x)
    vgy_new = 0; % no tangential velocity
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    %***Find x_new, y_new (note: a_new(3) is theta)
    X_new = X_curr + cos(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) - sin(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*timePerIter + 0.5*a_new(1)*timePerIter^2) + cos(a_new(3))*(vgy_new*timePerIter + 0.5*a_new(2)*timePerIter^2);
    X_curr = X_new;
    Y_curr = Y_new;
    
    %***Find change in theta between previous iterations (used to improve 
    %accuracy of finding t_curr)
    theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_change = 0;
    theta_prev = a_new(3);
    
    X_contact = X_curr - ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %For explination, see segment 1 description of this part
    
 
    %Find t_curr using sx_ge(6) equation (simplified) and X_contact
    t_curr = X_contact - 108/125 
    
    
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(7),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter

results_low = importdata('results_all_lo2.mat');
results_high = importdata('results_all_hi2.mat');

%find linear acceleration

%plot all the things
%plot xy cg points
hold on
subplot(3,3,[1 2 4 5]);
plot(results(1,:),results(2,:));
scatter(results(1,:),results(2,:));
title('Ball Position');

%Plot velocity over iterations
subplot(3,3,3)
hold on
plot(results(3,:))
plot(results_low(3,:))
plot(results_high(3,:))
title('Linear Velocity');

%plot acceleration
subplot(3,3,6)
hold on
plot(results(4,:))
plot(results_low(4,:))
plot(results_high(4,:))
title('Linear Acceleration');

%plot angular velocity
subplot(3,3,7)
hold on
plot(results(5,:))
plot(results_low(5,:))
plot(results_high(5,:))
title('Angular Velocity');

%plot angular acceleration
subplot(3,3,8)
hold on
plot(results(6,:))
plot(results_low(6,:))
plot(results_high(6,:))
title('Angular Acceleration');

%##########
%save('results_all_lo2.mat','results')
%##########



subplot(3,3,9)
hold on
plot(results(3,:))
plot(results_high(3,:))
plot(results_low(3,:))
title('thing');

%x = subs(evalin('base','sx'),t,t_cur); % eval. x at ti

% following code doesn't work atm, may try again later
% %y = subs(evalin('base','sy'),t,t_cur); % eval. y at ti
% 
% %Finds t_curr using the sx_ge equation if theta < 45, but switches to y
% %if > 45 deg for accuracy. t_curr +0.001 to prevent div. by zero errors
% 
% if abs(eval(subs(theta_ge(i),t,t_curr+0.001))) < pi/4  %if absolute value of theta is > than 45 deg;
%     t_curr = eval(solve(-X_curr + (75*tvar)/2 - (75*sin(tvar))/2 + 179/10)) %uses sx_ge(1) equation
% else
%     t_curr = eval(solve(-Y_curr + 75*cos(tvar)/2 + 844,'PrincipalValue',1)) % uses the sy_ge(1) equn
% end





%***** below is code used for freefall: unused in the final build?


% ***** Begin brief freefall of ball to segment 3
% for here, the only accel is due to gravity g
% omega remains constant and alpha = 0
% 
% %convert from funky normal & tangential velocity to plain ol' orthogonal
% vgx_new = cos(a_new(3))*curr_velo(1) - sin(a_new(3))*curr_velo(2);
% vgy_new = sin(a_new(3))*curr_velo(1) + cos(a_new(3))*curr_velo(2);
% curr_velo(1) = vgx_new;
% curr_velo(2) = vgy_new;
% 
% 
% 
% %Find (approximate) point of impact for when to switch to segment 3
% %this is done by guessing location of impact point via the figure
% t_curr = eval(solve(-0.2521+ (17*sin(t))/400 - (17*t)/400 + 51/200)); %(impact point x-coord approximated to 0.2521)
% X_impact = eval(subs(sx_ge(3),t,t_curr)) - sin(eval(subs(theta_ge(3),t,t_curr)))*ball_dim(2);
% %Approx. x-coords of ball center of mass when it impacts segment 3
% 
% while X_curr < X_impact
%     
%     %***Find v_new, omega, and angular accel
%     vgy_new = curr_velo(2) + -g*timePerIter; % accel down due to gravity
%     curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
%     %omega remains unchanged
%     curr_velo(4) = 0; % alpha = 0, no moment on ball
%     
%     %***Find x and y
%     X_curr = X_curr + curr_velo(1)*timePerIter;
%     Y_curr = Y_curr + curr_velo(2)*timePerIter;
%     results(1,iter) = X_curr; % x position using global coords
%     results(2,iter) = Y_curr; % y position using global coords
%     iter = iter + 1;
% end %end of free fall
%convert back from plain oatmeal orthogonal velocities to spicy sriracha
% %normal and tangential velocities
% vgx_new = cos(eval(subs(theta_ge(3),t,t_curr)))*curr_velo(1) + sin(eval(subs(theta_ge(3),t,t_curr)))*curr_velo(2);
% %vgy_new = - sin(eval(subs(theta_ge(3),t,t_curr)))*curr_velo(1) + cos(eval(subs(theta_ge(3),t,t_curr)))*curr_velo(2);
% vgy_new = 0; %***assume the ball doesn't bounce, and instead 'sticks' to the surface
% curr_velo(1) = vgx_new;
% curr_velo(2) = vgy_new;
% X_curr = X_impact; %'snap' to position on track
