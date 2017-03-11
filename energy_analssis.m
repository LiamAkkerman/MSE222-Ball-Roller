%TODO add into info

syms h_i h_f
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
h_i(i) = subs(sy(i), t, tmin(i));
h_f(i) = subs(sy(i), t, tmax(i));
%s(i) = int(d_arc_length(i), t, tmin(i), tmax(i));
%v(i);

%TODO make nomal force equation
work_friction(i) = int(normal_F(i)*d_arc_length(i), t, tmin(i), tmax(i));
work_gravity(i) = -int(ball_mass*g*(sy(i) - h_i(i)), t, tmin(1), tmax(i));

ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(omega_p*ball_dia/2)^2;


