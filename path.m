clear all;
clf;
%TODO add info


%basic mathmatical symbols
syms x y t h_brach s
%parametric equation symbols, to be array
syms sx sy  
%equations of Brachistochrone curve
syms sx_brach sy_brach  
syms arc_length

%project constants
hmax = 0.889; %starting hieght in m
%EDIT adding another segment so i form 6 to 7
segs = 7;   %number of line segments
m_in_inch = 0.0254; %how many m in an inch
ball_dia = 0.0081; %diameter of marble, in metres

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
rectangle('Position', [0 0 0.9144 0.9144]);
hold on;
circ_start_pos = [(m_in_inch - ball_dia) (hmax - ball_dia) 2*ball_dia 2*ball_dia];
circ_end_pos = [(hmax - ball_dia) (m_in_inch - ball_dia) 2*ball_dia 2*ball_dia];
rectangle('Position', circ_start_pos, 'Curvature',[1 1]);
rectangle('Position', circ_end_pos, 'Curvature',[1 1]);

i = 1; %let's try this way, still need to iterate
%substitute a values into brach functions
%the max is made up to have a upwards slope for the sick jump
tmin(i) = 0;
tmax(i) = 3.14*1.25;
h(i) = 0.075;
sx(i) = subs(sx_brach, h_brach, h(1)) + m_in_inch - ball_dia;
sy(i) = subs(sy_brach, h_brach, h(1)) + hmax - ball_dia;

%next line needs this
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%linear segment, maybe we'll keep, maybe not
i = 2;
%I made the max up randomly
%EDIT extend ramp so visable
tmin(i) = 0;
tmax(i) = 0.020;

%parametric equation of line
sx(i) = t + subs(sx(i-1), t, tmax(i-1));
sy(i) = t*subs(deriv(i-1), t, tmax(i-1)) + subs(sy(i-1), t, tmax(i-1));

%add it to the plot
ezplot(sx(i),sy(i),[tmin(i),tmax(i)])

%EDIT add aditional overhang to get rid of free fall between sick jump and
 %next segment
 i = 3;
 %find starting angel based on slope of sick jump
 tmin(i)=(pi/2)-atan(subs(deriv(i-2),t,tmax(i-2)));
 tmax(i)=pi;
 
 %Slop to bring ball form sick jump to second fast ramp
 h(i) = 2.5*ball_dia;
 sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1)) - (ball_dia)*cos(atan(subs(deriv(i-2),t,tmax(i-2)))+(pi/2));
 sy(i) = h(i)*sin(t) + subs(sy(i-1), t, tmax(i-1)) + (ball_dia)*cos(atan(subs(deriv(i-2),t,tmax(i-2)))+(pi/2));
%plot
ezplot(sx(i),sy(i),[tmin(i),tmax(i)]);

%EDIT cange all i to i+1 since adding another segment as the new i=3
%landing ramp
i = 4;
h(i) = 0.085;
tmin(i) = 0;
tmax(i) = 3.14;
%translational offset is purely guessed. we could caluculate better based on
%tridectory anaylsis
sx(i) = -subs(sx_brach, h_brach, h(i)) + subs(sx(i-1),t,tmax(i-1));
sy(i) = subs(sy_brach, h_brach, h(i)) + subs(sy(i-1),t,tmax(i-1));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%arc to bring the ball down to the last curve
i = 5;
h(i) = 5*ball_dia;
tmin(i) = 0;
tmax(i) = 3.14/2;
sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1));
sy(i) = h(i)*sin(t) + subs(sy(i-1 ), t, tmax(i-1)) - (h(i) - 2*ball_dia) + 0.001; %additional 1mm is for tolerance

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


%last large brach curve
i = 6;
h(i) = subs(sy(i-1), t, tmin(i-1)) - circ_end_pos(2);
tmin(i) = 0;
tmax(i) = 3.14;
sx(i) = subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));
sy(i) = subs(sy_brach, h_brach, h(i)) + subs(sy(i-1), t, tmin(i-1));

%standard brach curve is too long, scale it back.
%possibley change to scale both x and y instead of just x
x_factor = (circ_end_pos(1) + ball_dia - subs(sx(i-1), t, tmin(i-1)) - 0.025)/(subs(sx(i), t, tmax(i)) - subs(sx(i-1), t, tmin(i-1)));
sx(i) = x_factor*subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 
                        

%last line, flat
i = 7;
tmin(i) = 0;
tmax(i) = circ_end_pos(1) - subs(sx(i-1), t, tmax(i-1)) + ball_dia + m_in_inch;
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
    d_arc_length(i) = sqrt(diff(sy(i))^2 + diff(sx(i))^2)
    %following line depends on tmin = 0
    arc_length(i) = int(d_arc_length(i),t)
    s_arc_length(i) = solve(s == arc_length(i),t)
    %currently giving error. should output https://www.wolframalpha.com/input/?i=s%3D-(3+sin(t))%2F(40+sqrt(sin%5E2(t%2F2))),+solve+for+t
end