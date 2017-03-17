%general_test

%TODO
%Iterate over time via small step size & constant accel approximation.
%Double check if all units are right, or if we should convert to metres,
%not cm?

%Can fix approximation error by snapping to correct y-position given
%x-position if abs(theta) <= 45deg ? snap to x from y-pos if >45deg?

%Setup initial conditions
curr_velo = [0,0,0]; %vgx, vgy, omega
ball_dim = [0.000714,0.0081,1.8738216*10^-8]; %mass, radius, MoI
step = 0.1;  %step size (in seconds)

%Bring in info from path.m
tmin_ge = evalin('base','tmin'); %Initial t (vector) for each segment
tmax_ge = evalin('base','tmax'); %final t (vector) for each segment
t_curr = tmin_ge(1); %Starts at segment 1, so find initial value for 1st seg

sx_ge = evalin('base','sx'); %x-parametric equation
sy_ge = evalin('base','sy'); %x-parametric equation
X_curr = eval(subs(sx_ge,t,t_curr));
X_curr = X_curr(1); %Finds initial x,y coords
Y_curr = eval(subs(sy_ge,t,t_curr));
Y_curr = Y_curr(1); %Probably an easier way to do this lol



% NOTE: lowercase x,y coordinates represent the tangential and normal
% components of the ball for velocity, acceleration, etc analysis
% UPPERCASE X,Y coordinates represent the standard X Y coords (X = left to
% right, Y = down to up) for position analysis!



while t_curr < tmax_ge(1)
    %Find agx, agy:
    a_new = general_analysis(curr_velo,ball_dim, 1,0.5); %finds agx,agy, and 
    %theta given velocity, ball dimensions, segment #, and current t value
    
    %Find v_new, omega
    vgx_new = curr_velo(1) + a_new(1)*step;
    vgy_new = curr_velo(2) + a_new(2)*step;
    curr_velo(1) = vgx_new;
    curr_velo(2) = vgy_new; %this can be simplified to just reassign a variable to itself again?
    %TODO: add omega calcs
    
    %Find x_new, y_new
    X_new = X_curr + cos(a_new(3))*(vgx_new*step + a_new(1)*step^2) - sin(a_new(3))*(vgy_new*step + a_new(2)*step^2);
    Y_new = Y_curr + sin(a_new(3))*(vgx_new*step + a_new(1)*step^2) + cos(a_new(3))*(vgy_new*step + a_new(2)*step^2);
    X_curr = X_new;
    Y_curr = Y_new;
    %This is done by converting vgx & vgy tangential and normal velocities
    %(they should probably called something else lol) to the horizontal &
    %vertical coord system to find X_new and Y_new
    
    %Find new 't_curr'
    t_curr = solve(-X_curr + (75*tvar)/2 - (75*sin(tvar))/2 + 179/10) %uses sx_ga(1) equation
    
end


%x = subs(evalin('base','sx'),t,t_cur); % eval. x at ti
%y = subs(evalin('base','sy'),t,t_cur); % eval. y at ti