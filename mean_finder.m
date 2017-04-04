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
timevar = [0:0.002:0.948]

%plot all the things
%plot xy cg points
subplot(3,3,[1 2 4 5]);
hold on
plot(mean(1,:),mean(2,:),'r','LineWidth',1);
scatter(mean(1,scattervar),mean(2,scattervar),'k');
title('Ball Position');

%Plot velocity over iterations
subplot(3,3,3)
hold on
plot(mean(3,:),'LineWidth',1)
plot(results_low(3,:),'LineWidth',1)
plot(results_high(3,:),'LineWidth',1)
title('Linear Velocity');
legend('Mean','Min','Max')

%plot acceleration
subplot(3,3,6)
hold on
plot(mean(4,:),'LineWidth',1)
plot(results_low(4,:),'LineWidth',1)
plot(results_high(4,:),'LineWidth',1)
title('Linear Acceleration');
legend('Mean','Min','Max')

%plot angular velocity
subplot(3,3,7)
hold on
plot(mean(5,:),'LineWidth',1)
plot(results_low(5,:),'LineWidth',1)
plot(results_high(5,:),'LineWidth',1)
title('Angular Velocity');
legend('Mean','Min','Max')

%plot angular acceleration
subplot(3,3,8)
hold on
plot(mean(6,:),'LineWidth',1)
plot(results_low(6,:),'LineWidth',1)
plot(results_high(6,:),'LineWidth',1)
title('Angular Acceleration');
legend('Mean','Min','Max')
