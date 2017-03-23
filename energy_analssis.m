%uses symbolic euqations to calculate velocity at any given position
%currently in terms of t of each segment

syms y h
syms work_friction work_gravity normal_F ke v_g omega_p
%normal_F = nomal force a given posiion
%ke = kennetic enegy
%v_g = vlocity at point g (centre of mass)
%omega_p = angular veolicty where ball touches the track
%       is this right?

work_friction = sym('work_friction', [1,segs]);
work_gravity = sym('work_gravity', [1,segs]);
work_done = sym('work_done', [1,segs]);
ke = sym('ke', [1,segs]);
v_g = sym('v_g', [1,segs]);
h = sym('h', [1,segs]);
h_initial = NaN([1,segs]);

ball_mass = 0.000714; %in kgs
I_g = (2/5)*ball_mass*ball_radius^2;
g = 9.81;


%iterate to calculate heights for each segement
for i = 1:segs
    %change for if ball is on the underside of curve
    switch i
        case {1,2,3,5,6}
            %equation for hieght at any given point on curve
            h(i) = sy(i) + ball_radius*cos(theta(i));
            %initial hieght, used limit to avoid div by 0 errors
            h_initial(i) = limit(h(i), t, tmin(i));
        case {4}
            h(i) = sy(i) - ball_radius*cos(theta(i));
            h_initial(i) = limit(h(i), t, tmin(i));
    end
end

%iterate to calculate work for each segement
for i = 1:segs
    %neglecting firction currently
    work_friction(i) = 0;
    %integrates work done by gravity from the start of the curve to an abitrary t
    work_gravity(i) = -int(ball_mass*g, y, h_initial(i), h(i));
    %sums all work done
    work_done(i) = work_gravity(i) + work_friction(i);
end


%iterate to calculate kinneic energies for each segement
for i = 1:segs
    %change for rotation direction
    switch i
        %intial case with no intial KE. if vector index could be 0, his case would be unnesicary
        case {1}
            %kinnetic energy equation
            ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
            %solves energy equation to terms of t. quadradic so results in two values
            v_temp_vector = solve(0 + work_done(i) == ke(i), v_g(i));          
        %regular curve, postaive k rotation
        case {2,4,5,6}
            ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
            v_temp_vector = solve(subs(ke(i-1),t,tmax(i-1)) + work_done(i) == ke(i), v_g(i));
        %curves with negative k rotation
        case {3}
            %TODO figure out negative rotation
            ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
            v_temp_vector = solve(subs(ke(i-1),t,tmax(i-1)) + work_done(i) == ke(i), v_g(i));
    end
    
    %only take one result, and make it posative
    v_g(i) = abs(v_temp_vector(1));
    %reassign values of kinnetic energy now that v_g is known
    ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
end



i = 1;
%TODO not neglect friction
%normal_F(i) = ball_mass*g - (v_g(i)^2)/(rad_of_curv(i) - (ball_dia/2))
%work_friction(i) = int(normal_F(i)*d_arc_length(i), t, tmin(i), tmax(i))

%kinnetic energy equation
% ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
% %v_temp_vector = solve(subs(ke(i-1),t,tmax(i-1)) + work_done(i) == ke(i), v_g(i));
% v_temp_vector = solve(0 + work_done(i) == ke(i), v_g(i));
% %only take one result, and make it posative
% v_g(i) = abs(v_temp_vector(1));
% %reassign values of kinnetic energy now that v_g is known
% ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;

%display velocity at the end of the segment
double(subs(v_g(i), t, tmax(i)))




i = 2;
% %kinnnetic energy sum
% ke(i) = (1/2)*ball_mass*v_g(i)^2 + (1/2)*I_g*(v_g(i)/ball_radius)^2;
% %solves energy equation to terms of t. quadradic so results in two values
% v_temp_vector = solve(subs(ke(i-1),t,tmax(i-1)) + work_done(i) == ke(i), v_g(i));
% %only take one result, and make it posative
% v_g(i) = abs(v_temp_vector(1));

%display velocity at the end of the segment
double(subs(v_g(i), t, tmax(i)))

