%TODO add into info

syms h_i h_f y
syms work_friction work_gravity normal_F ke v_g omega_p
%normal_F = nomal force a given posiion
%ke = kennetic enegy
%v_g = vlocity at point g (centre of mass)
%omega_p = angular veolicty where ball touches the track
%       is this right?

ball_mass = 0.0015; %in kgs
I_g = (2/5)*ball_mass*(ball_dia/2)^2;
g = 9.81;

i = 1;
%h_i(i) = subs(sy(i), t, tmin(i));
%h_f(i) = subs(sy(i), t, tmax(i));
%s(i) = int(d_arc_length(i), t, tmin(i), tmax(i));
%v(i);

%TODO make nomal force equation
%normal_F(i) = ball_mass*g - (v_g(i)^2)/(rad_of_curv(i) - (ball_dia/2))
%work_friction(i) = int(normal_F(i)*d_arc_length(i), t, tmin(i), tmax(i))
work_friction(i) = 0;
work_gravity(i) = -int(ball_mass*g, y, subs(sy(i), t, tmin(i)), sy(i))
work_done(i) =  work_gravity(i) + work_friction(i);

ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g/ball_dia/2)^2;
%solve(subs(ke(i-1),t,tmax(i-1)) + (work_fiction + work_graviy) == ke(i), v_g(i)); 
v_temp_vector = solve(0 + work_done(i) == ke(i), v_g(i))
v_g(i) = abs(v_temp_vector(1))

double(subs(v_g(i), t, tmax(i)))


i = 2;
sy(i)
work_friction(i) = 0;
work_gravity(i) = -int(ball_mass*g, y, subs(sy(i), t, tmin(i)), sy(i))
work_done(i) =  work_gravity(i) + work_friction(i)

ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g/ball_dia/2)^2
v_temp_vector = solve(subs(ke(i-1),t,tmax(i-1)) + work_done(i) == ke(i), v_g(i)) 
%v_temp_vector = solve(0 + work_done == ke(i), v_g(i));
v_g(i) = abs(v_temp_vector(1))

double(subs(v_g(i), t, tmax(i)))

