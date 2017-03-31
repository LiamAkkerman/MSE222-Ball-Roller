function [ results ] = recordresults( X_curr, Y_curr, curr_velo, a_new, iter )
% Records results to an array (cleaning up the main code)

    results(1) = X_curr; % x position using global coords
    results(2) = Y_curr; % y position using global coords
    results(3) = abs(curr_velo(1)); % magnitude of velocity (is actually just tangential velo, as no normal component of velocity)
    results(4) = sqrt(a_new(1)^2 + a_new(2)^2); % magnitude of accel
    results(5) = curr_velo(3); % angular velocity
    results(6) = curr_velo(4); % angular acceleration
    results(7) = a_new(3); %theta
    results(8) = a_new(1); %tangential accel
    results(9) = a_new(2); %normal accel
    %results(10) = current friction force?




end

