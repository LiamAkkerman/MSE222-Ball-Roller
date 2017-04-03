function [ soln ] = kinn_anal(anal_param,ball_dim,path_array,theta,rad_of_curv)

%Setup initial conditions
curr_velo = [0,0,0,0]; %vgx, vgy, omega, alpha
%ball_dim = [0.000551,0.0081,1.8738216*10^-8]; %mass (kg), radius (m), MoI (kg/m^2)
timePerIter = anal_param(1); %0.005;  %step size (in seconds)
atan_tol = anal_param(2); %the tolerance given to atan when limit does not exist
iter = 1;
syms t tvar;
g = 9.81;

%Bring in info from path.m
tmin = eval(path_array(3,:)); %Initial t (vector) for each segment
tmax = eval(path_array(4,:)); %final t (vector) for each segment
t_curr = tmin(1); %Starts at segment 1, so find initial value for 1st seg

sx = path_array(1,:); %x-parametric equation
sy = path_array(2,:); %y-parametric equation

X_curr = eval(subs(sx,t,t_curr));
X_curr = X_curr(1) + ball_dim(2)*cos(eval(subs(theta(1),t,t_curr+0.001))+pi/2); %Finds initial x,y coords of the center of mass
Y_curr = eval(subs(sy,t,t_curr));
Y_curr = Y_curr(1) + ball_dim(2)*sin(eval(subs(theta(1),t,t_curr+0.001))+pi/2);
%this is found by finding the position of contact on the curve, then adding
%the radius of the ball to find the CoM (radius is in m, so convert to mm)
theta_prev = eval(subs(theta(1),t,t_curr+0.001)); %used later for accuracy improvement for finding X_contact




% NOTE: lowercase x,y coordinates represent the tangential and normal
% components of the ball for velocity, acceleration, etc analysis
% UPPERCASE X,Y coordinates represent the global X Y coords (X = left to
% right, Y = down to up) for position analysis!

% ***** First segment
while t_curr < tmax(1)
    disp('segment 1');
    %***Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 1,t_curr, atan_tol); %finds agx,agy, and 
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
    Y_curr = eval(subs(sy(1),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
   
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
t_curr = tmin(2); %this isn't always true for the start of each segment!

while t_curr < tmax(2)
    disp('segment 2');
    %***Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 2,t_curr, atan_tol); %finds agx,agy, and 
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
    t_curr = eval(solve(-X_contact + tvar - (3*sin(157/40))/80 + 13159/80000)) %uses sx_ge(2) equation
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy(2),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
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
t_curr = tmin(3)

while t_curr < tmax(3)
    disp('segment 3');

    % *** Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 3,t_curr, atan_tol); %finds agx,agy,theta
    
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
    t_curr = acos(102823038966022130827/9727775195120271360 - (4000*X_contact)/81)
    %Above equation was found by solving the sx_ge(3) == X_contact equation for t 
    
    % *** Snap Y_curr to line to account for slight drift 
    Y_curr = eval(subs(sy(3),t,t_curr)) + ball_dim(2)*sin(a_new(3)-pi/2);
   
    % *** Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1; 
end
time_taken = (iter-1)*timePerIter


% ***** Fourth Segment
t_curr = tmin(4); %this isn't always true for the start of each segment!

curr_velo(1) = -curr_velo(1); % To account for ball travelling right to left
theta_prev = 1.5707; % To account for theta_change being unused in previous segment

while t_curr < tmax(4)
    disp('segment 4');

    %***Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 4,t_curr, atan_tol); %finds agx,agy, and 
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
    Y_curr = eval(subs(sy(4),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter


% ***** Segment 5

t_curr = tmax(5)

%NOTE!!: For this segment, t goes from tmax to tmin for some reason!

while t_curr > tmin(5)
    disp('segment 5');

    % *** Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 5,t_curr, atan_tol); %finds agx,agy,theta
    
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
    %Above equation was found by solving the sx_ge(3) == X_contact equation for t 
    
    % *** Snap Y_curr to line to account for slight drift 
    Y_curr = eval(subs(sy(5),t,t_curr)) + ball_dim(2)*sin(a_new(3)-pi/2);
   
    % *** Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1; 
end



% ***** Segment 6
t_curr = tmin(6); %this isn't always true for the start of each segment!

curr_velo(1) = -curr_velo(1); % To account for ball travelling left to right
theta_prev = -1.5704;

while t_curr < tmax(6)
    disp('segment 6');

    %***Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 6,t_curr, atan_tol); %finds agx,agy, and 
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
    t_curr = solve((2520299339554973*t)/9843000955906036 - (2520299339554973*sin(t))/9843000955906036 + 2176540528804513/36028797018963968 - X_contact)
    
    
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy(6),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter

% ***** Segment 7
t_curr = tmin(7); %this isn't always true for the start of each segment!

while X_curr < 0.889 % While ball position isn't at the finish position
    disp('segment 7');

    %***Find agx, agy:
    a_new = general_analysis(theta,rad_of_curv,curr_velo,ball_dim, 7,t_curr, atan_tol); %finds agx,agy, and 
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
    Y_curr = eval(subs(sy(7),t,t_curr)) + ball_dim(2)*sin(a_new(3)+pi/2);
    
    %***Save values in array
    results(:,iter) = recordresults(X_curr,Y_curr,curr_velo,a_new,iter);
    
    
    
    iter = iter + 1;
end
time_taken = (iter-1)*timePerIter



%plot all the things
%plot xy cg points
hold on
subplot(3,3,[1 2 4 5]);
plot(results(1,:),results(2,:));
scatter(results(1,:),results(2,:));
title('Ball Position');

%Plot velocity over iterations
subplot(3,3,3)
plot(results(3,:))
title('Linear Velocity');

%plot acceleration
subplot(3,3,6)
plot(results(4,:))
title('Linear Acceleration');

%plot angular velocity
subplot(3,3,7)
plot(results(5,:))
title('Angular Velocity');

%plot angular acceleration
subplot(3,3,8)
plot(results(6,:))
title('Angular Acceleration');

end

