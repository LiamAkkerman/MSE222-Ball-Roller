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
segs = 7;   %number of line segments
mmininch = 25.4; %how many mm in an inch
balldia = 15; %diameter of marble, 15 is a guess

%initialize empty vectors for hieghts, starts, and ends
h = NaN(segs);
tmin = NaN(segs);
tmax = NaN(segs);

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
circ_start_pos = [(mmininch - 0.5*balldia) (hmax - 0.5*balldia) balldia balldia];
circ_end_pos = [(hmax - 0.5*balldia) (mmininch - 0.5*balldia) balldia balldia];
rectangle('Position', circ_start_pos, 'Curvature',[1 1]);
rectangle('Position', circ_end_pos, 'Curvature',[1 1]);

i = 1; %let's try this way, still need to iterate
%substitute a values into brach functions
h(i) = 75;
sx(i) = subs(sx_brach, h_brach, h(1)) + mmininch - 0.5*balldia;
sy(i) = subs(sy_brach, h_brach, h(1)) + hmax - 0.5*balldia;
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));
%TODO make function for all of these

%the max is made up to have a upwards slope for the sick jump
tmin(i) = 0;
tmax(i) = 3.14*1.25;
%TODO automate max t

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 

i = 2;
%linear segment, maybe we'll keep, maybe not
%inherits perovous derivitive
deriv(i) = subs(deriv(i-1), t, tmax(i-1));
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));

%parametric equation of line
sx(i) = t + subs(sx(i-1), t, tmax(i-1));
sy(i) = t*deriv(i) + subs(sy(i-1), t, tmax(i-1));

%I made the max up randomly
tmin(i) = 0;
tmax(i) = 10;

%add it to the plot
ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


i = 3;
h(i) = 85;
sx(i) = -subs(sx_brach, h_brach, h(i)) + 255;
sy(i) = subs(sy_brach, h_brach, h(i)) + 839;
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));

tmin(i) = 0;
tmax(i) = 3.14;

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


i = 4;
h(i) = 3*balldia;
tmin(i) = 0;
tmax(i) = 3.14/2;
sx(i) = -h(i)*cos(t) + subs(sx(i-1), t, tmax(i-1));
sy(i) = h(i)*sin(t) + subs(sy(i-1), t, tmax(i-1)) + 1; %additional 1 is for tolerance

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 
                        
%for i = 1:segs
%    ezplot(sx(i),sy(i),[tmin(i),tmax(i)])     %plot sections   
%end    

