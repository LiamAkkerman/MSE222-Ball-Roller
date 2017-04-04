m1 = importdata('results_mass_hi2.mat');
m2 = importdata('results_mass_lo2.mat');
r1 = importdata('results_radius_hi2.mat');
r2 = importdata('results_radius_lo2.mat');
f1 = importdata('results_frict_hi2.mat');
f2 = importdata('results_frict_lo2.mat');

results = importdata('results_standard2.mat');
results_high = importdata('results_all_hi2.mat');
results_low  = importdata('results_all_lo2.mat');


for j = 1:475
    for i = 1:9
        mean(i,j) = (m1(i,j) + m2(i,j) + r1(i,j) + r2(i,j) + f1(i,j) + f2(i,j))/6;
    end
end
%plot(mean(1,:),mean(2,:))

scattervar = [1:10:length(mean)];

%plot all the things
%plot xy cg points
subplot(3,3,[1 2 4 5]);
hold on
plot(mean(1,:),mean(2,:),'r');
scatter(mean(1,scattervar),mean(2,scattervar),'k');
title('Ball Position');

%Plot velocity over iterations
subplot(3,3,3)
hold on
plot(mean(3,:))
plot(results_low(3,:))
plot(results_high(3,:))
title('Linear Velocity');

%plot acceleration
subplot(3,3,6)
hold on
plot(mean(4,:))
plot(results_low(4,:))
plot(results_high(4,:))
title('Linear Acceleration');

%plot angular velocity
subplot(3,3,7)
hold on
plot(mean(5,:))
plot(results_low(5,:))
plot(results_high(5,:))
title('Angular Velocity');

%plot angular acceleration
subplot(3,3,8)
hold on
plot(mean(6,:))
plot(results_low(6,:))
plot(results_high(6,:))
title('Angular Acceleration');
