function [ soln ] = rad_of_curv_funct(sx,sy)
    size = length(sx);
    soln = sym('rad_of_curv', [1,size]);
    for i = 1:size
        soln(i) = simplify(((diff(sx(i))^2 + diff(sy(i))^2)^(3/2)) / abs( (diff(sx(i))*(diff(sy(i),2))) - (diff(sy(i))*(diff(sx(i),2))) )); 
    end
end


