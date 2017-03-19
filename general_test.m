%general_test

%TODO
% add code for each of the segments
% NOTE: if any of the original path.m equations are changed, the code for
% finding the next t_curr must be fixed as well as it is hardcoded in atm

%Can fix approximation error by snapping to correct y-position given
%x-position if abs(theta) <= 45deg ? snap to x from y-pos if >45deg?

%Setup initial conditions
curr_velo = [0,0,0,0]; %vgx, vgy, omega, alpha
ball_dim = [0.000714,0.0081,1.8738216*10^-8]; %mass, radius, MoI
step = 0.01;  %step size (in seconds)
iter = 1;
syms tvar;

%Bring in info from path.m
tmin_ge = evalin('base','tmin'); %Initial t (vector) for each segment
tmax_ge = evalin('base','tmax'); %final t (vector) for each segment
t_curr = tmin_ge(1); %Starts at segment 1, so find initial value for 1st seg

sx_ge = evalin('base','sx'); %x-parametric equation
sy_ge = evalin('base','sy'); %x-parametric equation
theta_ge = evalin('base','theta');

X_curr = eval(subs(sx_ge,t,t_curr));
X_curr = X_curr(1) +1000*ball_dim(2)*cos(eval(subs(theta_ge(1),t,t_curr+0.001))+pi/2); %Finds initial x,y coords of the center of mass
Y_curr = eval(subs(sy_ge,t,t_curr));
Y_curr = Y_curr(1)+1000*ball_dim(2)*sin(eval(subs(theta_ge(1),t,t_curr+0.001))+pi/2);
%this is found by finding the position of contact on the curve, then adding
%the radius of the ball to find the CoM (radius is in m, so convert to mm)
theta_prev = eval(subs(theta_ge(1),t,t_curr+0.001)); %used later for accuracy improvement




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
    vgx_new = curr_velo(1) + a_new(1)*step; %tangential (x)
    vgy_new = curr_velo(2) + a_new(2)*step; %normal (y)
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    
    %***Find X_new, Y_new
    X_new = 1000*(X_curr/1000 + cos(a_new(3))*(vgx_new*step + a_new(1)*step^2) - sin(a_new(3))*(vgy_new*step + a_new(2)*step^2));
    Y_new = 1000*(Y_curr/1000 + sin(a_new(3))*(vgx_new*step + a_new(1)*step^2) + cos(a_new(3))*(vgy_new*step + a_new(2)*step^2));
    X_curr = X_new;
    Y_curr = Y_new; %redundant, can be removed?
    %This is done by converting vgx & vgy tangential and normal velocities
    %(they should probably called something else lol) to the horizontal &
    %vertical coord system to find X_new and Y_new. X & Y are in mm, so
    %convert to meters first.
    
    %***Find change in theta between previous iterations (used to improve 
    %accuracy of finding t_curr)
    theta_change = a_new(3) - theta_prev; %(i-1 theta) - (i-2 theta) to approx change in theta
    theta_prev = eval(subs(theta_ge(1),t,t_curr+0.001)); % +0.001 added to fix div by zero errors
    
    %***Find new 't_curr' by calculating the x-coord of the point of
    %contact
    
    X_contact = X_curr - 1000*ball_dim(2)*cos(a_new(3)+pi/2+theta_change); 
    %Some weird geometry stuff going on here. Basically, this accounts for
    %how the change of the position of the center of mass isn't exactly the
    %same as the change of the position of the point of contact. Pos. of
    %the contact point is used for finding the next t_curr value, so it's
    %fairly important to account for the difference.
    
    t_curr = eval(solve(-X_contact + (75*tvar)/2 - (75*sin(tvar))/2 + 179/10)) %uses sx_ge(1) equation
    % 'snap' Y_curr to curve + add radius of ball to account for iteration error
    Y_curr = eval(subs(sy_ge(1),t,t_curr)) + 1000*ball_dim(2)*sin(a_new(3)+pi/2);
   
    %***Save values in array
    results(1,iter) = X_curr; % x position using global coords
    results(2,iter) = Y_curr; % y position using global coords
    results(3,iter) = (curr_velo(1)^2+curr_velo(2)^2)^0.5; % magnitude of velo
    results(4,iter) = (a_new(1)^2 + a_new(2)^2)^0.5; % magnitude of accel
    results(5,iter) = curr_velo(3); % angular velocity
    results(6,iter) = curr_velo(4); % angular acceleration
    
    iter = iter + 1;
end
time_taken = (iter-1)*step


% ***** Second Segment
t_curr = tmin_ge(2); %this isn't always true for the start of each segment!

while t_curr < tmax_ge(2)
    %***Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 2,t_curr); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %***Find v_new, omega, and angular accel
    vgx_new = curr_velo(1) + a_new(1)*step; %tangential (lowercase x)
    vgy_new = curr_velo(2) + a_new(2)*step; %normal (lowercase y)
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    curr_velo(3) = vgx_new/ball_dim(2); % omega = v/r
    curr_velo(4) = a_new(1)/ball_dim(2); % alpha = agx/r
    
    
    %***Find x_new, y_new
    X_new = 1000*(X_curr/1000 + cos(a_new(3))*(vgx_new*step + a_new(1)*step^2) - sin(a_new(3))*(vgy_new*step + a_new(2)*step^2));
    Y_new = 1000*(Y_curr/1000 + sin(a_new(3))*(vgx_new*step + a_new(1)*step^2) + cos(a_new(3))*(vgy_new*step + a_new(2)*step^2));
    X_curr = X_new;
    Y_curr = Y_new;
    %This is done by converting vgx & vgy tangential and normal velocities
    %(they should probably called something else lol) to the horizontal &
    %vertical coord system to find X_new and Y_new. X & Y are in mm, so
    %convert to meters first.
    
    %***Find new 't_curr' using x coords
    
    t_curr = eval(solve(-X_curr + tvar - (75*sin(157/40))/2 + 13207/80)) %uses sx_ge(2) equation
    % 'snap' Y_curr to curve to account for iteration error
    Y_curr = eval(subs(sy_ge(2),t,t_curr));
    
    %***Save values in array
    results(1,iter) = X_curr; % x position using global coords
    results(2,iter) = Y_curr; % y position using global coords
    results(3,iter) = (curr_velo(1)^2+curr_velo(2)^2)^0.5; % magnitude of velo
    results(4,iter) = (a_new(1)^2 + a_new(2)^2)^0.5; % magnitude of accel
    results(5,iter) = curr_velo(3); % angular velocity
    results(6,iter) = curr_velo(4); % angular acceleration
    
    iter = iter + 1;
end
time_taken = (iter-1)*step

hold on
scatter(results(1,:),results(2,:))

% ***** Begin brief freefall of ball to segment 3
% for here, the only accel is due to gravity g
% omega remains constant and alpha = 0

% TODO add code here


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