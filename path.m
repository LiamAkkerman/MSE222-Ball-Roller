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



%define first section hieght
h(1) = 3*mmininch;        


i = 1; %let's try this way, still need to iterate
%substitute a values into brach functions
sx(i) = subs(sx_brach, h_brach, h(1));    %define parametic equations for 1st Brachistochrone curve
sy(i) = subs(sy_brach, h_brach, h(1)) + hmax;
deriv(i) = simplify(diff(sy(i))/diff(sx(i)));
deriv2(i) = simplify(diff(deriv(i)))
theta(i) = simplify(atan(deriv(i)))


tmin(i) = 0;
tmax(i) = 3.14*1.25;           %max t for 1st section
%TODO automate max t

double(subs(theta(i), t, tmax(i)))
                        

                        
%for i = 1:segs
    ezplot(sx(i),sy(i),[tmin(i),tmax(i)])     %plot sections   
%end    

