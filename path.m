clear all;
clf;
%TODO add info


%basic mathmatical symbols
syms x y t h_brach
%parametric equation symbols, to be array
syms sx sy  
%equations of Brachistochrone curve
syms sx_brach sy_brach  

%project constants
hmax = 889; %starting hieght in mm
segs = 6;   %number of line segments
mmininch = 25.4; %how many mm in an inch
ball_dia = 8.1; %diameter of marble, 15 is a guess

%initialize empty vectors for hieghts, starts, and ends
h = NaN([1,segs]);
tmin = NaN([1,segs]);
tmax = NaN([1,segs]);
%deriv = NaN(segs);
%deriv2 = NaN(segs);
%theta = NaN(segs);

%vectors of parameetric equations for each line segments
sx = sym('sx', [1,segs]);
sy = sym('sy', [1,segs]);
deriv = sym('sy', [1,segs]);

%define Brachistochrone curve parametric equations
sx_brach = 0.5*h_brach*(t-sin(t));  
sy_brach = -0.5*h_brach*(1-cos(t));


%plot borders
rectangle('Position', [0 0 914.4 914.4]);
hold on;
circ_start_pos = [(mmininch - 0.5*ball_dia) (hmax - 0.5*ball_dia) ball_dia ball_dia];
circ_end_pos = [(hmax - 0.5*ball_dia) (mmininch - 0.5*ball_dia) ball_dia ball_dia];
rectangle('Position', circ_start_pos, 'Curvature',[1 1]);
rectangle('Position', circ_end_pos, 'Curvature',[1 1]);

i = 1; %let's try this way, still need to iterate
%substitute a values into brach functions
%the max is made up to have a upwards slope for the sick jump
tmin(i) = 0;
tmax(i) = 3.14*1.25;
h(i) = 75;
sx(i) = subs(sx_brach, h_brach, h(1)) + mmininch - 0.5*ball_dia;
sy(i) = subs(sy_brach, h_brach, h(1)) + hmax - 0.5*ball_dia;

%next line needs this
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%linear segment, maybe we'll keep, maybe not
i = 2;
%I made the max up randomly
tmin(i) = 0;
tmax(i) = 10;

%parametric equation of line
sx(i) = t + subs(sx(i-1), t, tmax(i-1));
sy(i) = t*subs(deriv(i-1), t, tmax(i-1)) + subs(sy(i-1), t, tmax(i-1));

%add it to the plot
ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%landing ramp
i = 3;
h(i) = 85;
tmin(i) = 0;
tmax(i) = 3.14;
%translational offset is purely guessed. we could caluculate better based on
%tridectory anaylsis
sx(i) = -subs(sx_brach, h_brach, h(i)) + 255;
sy(i) = subs(sy_brach, h_brach, h(i)) + 839;

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%arc to bring the ball down to the last curve
i = 4;
h(i) = 2.5*ball_dia;
tmin(i) = 0;
tmax(i) = 3.14/2;
sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1));
sy(i) = h(i)*sin(t) + subs(sy(i-1 ), t, tmax(i-1)) - (h(i) - ball_dia) + 1; %additional 1 is for tolerance

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%last large brach curve
i = 5;
h(i) = subs(sy(i-1), t, tmin(i-1)) - circ_end_pos(2);
tmin(i) = 0;
tmax(i) = 3.14;
sx(i) = subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));
sy(i) = subs(sy_brach, h_brach, h(i)) + subs(sy(i-1), t, tmin(i-1));

%standard brach curve is too long, scale it back.
%possibley change to scale both x and y instead of just x
x_factor = (circ_end_pos(1) + 0.5*ball_dia - subs(sx(i-1), t, tmin(i-1)) - 25)/(subs(sx(i), t, tmax(i)) - subs(sx(i-1), t, tmin(i-1)));
sx(i) = x_factor*subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 
                        

%last line, flat
i = 6;
tmin(i) = 0;
tmax(i) = circ_end_pos(1) - subs(sx(i-1), t, tmax(i-1)) + 0.5*ball_dia + mmininch;
sx(i) = t + subs(sx(i-1), t, tmax(i-1));
sy(i) = subs(sy(i-1), t, tmax(i-1));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)])



%get more info for anaylis
for i = 1:segs
    deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
    deriv2(i) = simplify(diff(deriv(i)));
    theta(i) = simplify(atan(deriv(i)));
    
    %rad_of_curv(i) = ((1 + deriv(i)^2)^(3/2))/deriv2(i); wrong form of eq
    
    rad_of_curv(i) = simplify(((diff(sx(i))^2 + diff(sy(i))^2)^(3/2)) / abs( (diff(sx(i))*(diff(sy(i),2))) - (diff(sy(i))*(diff(sx(i),2))) )); 
    %need to use parametric equation for rad of curv (what a beastcyclo)
    d_arc_length(i) = simplify(sqrt(diff(sy(i))^2 + diff(sx(i))^2));
end