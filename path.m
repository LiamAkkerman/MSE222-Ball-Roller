function [ Soln ] = path(parameters,ball_dim)


start_hieght = parameters(1);
ramp_angle = parameters(2);
ramp_length = parameters(3);
arc_1_dia = parameters(4);
landing_height = parameters(5);
arc_2_dia = parameters(6);
flat_legnth = parameters(7);

%basic mathmatical symbols
syms x y t h_brach s
%parametric equation symbols, to be array
syms sx sy  
%equations of Brachistochrone curve
syms sx_brach sy_brach  
%syms arc_length

%project constants
hmax = 0.889; %starting hieght in m
segs = parameters(8);   %number of line segments
tolerance = parameters(9);   %tolerance between paths when ball switches sides
m_in_inch = 0.0254; %how many m in an inch
ball_dia = 2*ball_dim(2); %diameter of marble, in metres

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

% arc_length = sym('arc_length', [1,segs]);
% t_arc_length = sym('t_arc_length', [1,segs]);
% d_arc_length = sym('d_arc_length', [1,segs]);

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


%first segment, brach surve into ramp
i = 1;
    %the max is made up to have a upwards slope for the sick jump
    tmin(i) = 0;
    tmax(i) = pi*ramp_angle; %1.25
    h(i) = start_hieght; %0.075
    sx(i) = subs(sx_brach, h_brach, h(1)) + m_in_inch - ball_dia; %substitute a values into brach functions
    sy(i) = subs(sy_brach, h_brach, h(1)) + hmax - ball_dia;

    %next line needs this
    ramp_slope = subs((diff(sy(i))/diff(sx(i))), t, tmax(i));
    ramp_angle = atan(ramp_slope);



%linear segment that is the ramp
i = 2;
    tmin(i) = 0;
    tmax(i) = ramp_length; %0.020
    sx(i) = t + subs(sx(i-1), t, tmax(i-1));
    sy(i) = t*ramp_slope + subs(sy(i-1), t, tmax(i-1));


%arc to bring ball form sick jump to second fast ramp
i = 3;
    tmin(i)=(pi/2)-ramp_angle; %find starting angel based on slope of sick jump
    tmax(i)=pi;
    h(i) = arc_1_dia*ball_dia; %2.5
    sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1)) - (ball_dia)*cos(ramp_angle+(pi/2));
    sy(i) = h(i)*sin(t) + subs(sy(i-1), t, tmax(i-1)) + (ball_dia)*cos(ramp_angle+(pi/2));

 
%landing ramp
i = 4;
    h(i) = landing_height; %0.085
    tmin(i) = 0;
    tmax(i) = pi;
    sx(i) = -subs(sx_brach, h_brach, h(i)) + subs(sx(i-1),t,tmax(i-1));
    sy(i) = subs(sy_brach, h_brach, h(i)) + subs(sy(i-1),t,tmax(i-1));


%arc to bring the ball down to the last curve
i = 5;
    h(i) = arc_2_dia*ball_dia;
    tmin(i) = 0;
    tmax(i) = pi/2;
    sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1));
    sy(i) = h(i)*sin(t) + subs(sy(i-1 ), t, tmax(i-1)) - (h(i) - 2*ball_dia) + tolerance;


%last large brach curve
i = 6;
    h(i) = subs(sy(i-1), t, tmin(i-1)) - circ_end_pos(2);
    tmin(i) = 0;
    tmax(i) = 3.14;
    sx(i) = subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));
    sy(i) = subs(sy_brach, h_brach, h(i)) + subs(sy(i-1), t, tmin(i-1));

    %standard brach curve is too long, scale it back.
    %possibley change to scale both x and y instead of just x
    x_factor = (circ_end_pos(1) + ball_dia - subs(sx(i-1), t, tmin(i-1)) - flat_legnth)/(subs(sx(i), t, tmax(i)) - subs(sx(i-1), t, tmin(i-1)));
    sx(i) = x_factor*subs(sx_brach, h_brach, h(i)) + subs(sx(i-1), t, tmin(i-1));
                        

%last line, flat
i = 7;
    tmin(i) = 0;
    tmax(i) = circ_end_pos(1) - subs(sx(i-1), t, tmax(i-1)) + ball_dia + m_in_inch;
    sx(i) = t + subs(sx(i-1), t, tmax(i-1));
    sy(i) = subs(sy(i-1), t, tmax(i-1));


    
    
    
    
    
    
    

%get more info for anaylis
for i = 1:segs
    %following lines commented out for testing purposes
    deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
%     deriv2(i) = simplify(diff(deriv(i)));
    theta(i) = simplify(atan(deriv(i)));
    
    %rad_of_curv(i) = ((1 + deriv(i)^2)^(3/2))/deriv2(i); wrong form of eq
    
    %commented for speed
    rad_of_curv(i) = simplify(((diff(sx(i))^2 + diff(sy(i))^2)^(3/2)) / abs( (diff(sx(i))*(diff(sy(i),2))) - (diff(sy(i))*(diff(sx(i),2))) )); 
    %need to use parametric equation for rad of curv (what a beastcyclo)
%     d_arc_length(i) = sqrt(diff(sy(i))^2 + diff(sx(i))^2);
%     %following line depends on tmin = 0
%     if i == 6
%         digits(8);
%         d_arc_length(i) = vpa(d_arc_length(i));
%     end
%     if tmin(i) == 0
%         arc_length(i) = int(d_arc_length(i),t);
%     else
%         arc_length(i) = int(d_arc_length(i),t,tmin(i),t);
%     end
%     %t_arc_length2(i) = solve(s == arc_length(i),t,'ReturnConditions',1,'PrincipalValue',true)
%     gg = solve(s == arc_length(i),t,'Real',1)
%     %currently giving error. should output https://www.wolframalpha.com/input/?i=s%3D-(3+sin(t))%2F(40+sqrt(sin%5E2(t%2F2))),+solve+for+t
%     
%     string = ['segement = ', num2str(i), ',   ds = ', char(d_arc_length(i)), ',   s(t) = ', char(arc_length(i)), ',   t(s) = ', char(t_arc_length(i))];
%     disp(string);
end

end

%vpa(arc_length(6))