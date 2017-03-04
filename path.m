syms x y t
syms s1 s2 s3 s4 s5 s6 s7
syms s1x s1y            %parametric equation symbols

a1 = (3/2)*25.4;        %define first section hieght

s1x = a1*(t-sin(t));    %define parametic equations for 1st Brachistochrone curve
s1y = -a1*(1-cos(t));

%t1max = solve(s1y == -(3*25.4), t)
t1max = 3.14;           %max t for 1st section
                        %TODO automate max t

fplot(s1x,s1y+914.4-25.4,[0,t1max])     %plot sections
