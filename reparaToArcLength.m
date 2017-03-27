%uses parametric equations to reparametrize as arclength
syms Length; %arclength for a given segment
syms fun; %placeholder so intagration equation is readable
%vectors of arclength and fun for each segment
Length = sym('Length', [1,segs]);
fun = sym('fun', [1,segs]);

%start with first segment
i = 1;
%dx/dt is always positive so no changes needed
fun(i) = sqrt((diff(sx(i))^2 + diff(sy(i))^2))
Length(i) = integral(fun(i),tmin(i),tmax(i))
