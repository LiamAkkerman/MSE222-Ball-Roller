function [ soln ] = theta_funct(sx,sy)
    size = length(sx);
    soln = sym('theta', [1,size]);
    for i = 1:size
        soln(i) = atan(diff(sy(i))/diff(sx(i)));
    end
end
