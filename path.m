clear all;
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


     


i = 1; %let's try this way, still need to iterate
%substitute a values into brach functions
h(i) = 75;
sx(i) = subs(sx_brach, h_brach, h(1));    %define parametic equations for 1st Brachistochrone curve
sy(i) = subs(sy_brach, h_brach, h(1)) + hmax;
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));
%TODO make function for all of these

%the max is made up to have a upwards slope for the sick jump
tmin(i) = 0;
tmax(i) = 3.14*1.25;
%TODO automate max t

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 
hold on;

i = 2;
%linear segment, maybe we'll keep, maybe not
%inherits perovous derivitive
deriv(i) = subs(deriv(i-1), t, tmax(i-1));
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));

%parametric equation of line
sx(i) = subs(sx(i-1), t, tmax(i-1)) + t;
sy(i) = subs(sy(i-1), t, tmax(i-1)) + t*deriv(i);

%I made the max up randomly
tmin(i) = 0;
tmax(i) = 10;

%add it to the plot
ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 


i = 3;
h(i) = 85;
sx(i) = -subs(sx_brach, h_brach, h(i)) + 230;    %define parametic equations for 1st Brachistochrone curve
sy(i) = subs(sy_brach, h_brach, h(i)) + 850;
deriv(i) = 0.45;
deriv2(i) = simplify(diff(deriv(i)));
theta(i) = simplify(atan(deriv(i)));

tmin(i) = 0;
tmax(i) = 3.14;

ezplot(sx(i),sy(i),[tmin(i),tmax(i)]) 
                        
%for i = 1:segs
%    ezplot(sx(i),sy(i),[tmin(i),tmax(i)])     %plot sections   
%end    

