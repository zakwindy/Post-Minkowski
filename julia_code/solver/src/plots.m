clear;
close all;

data = importdata('PMdata.csv');

x1 = data(1,:);
y1 = data(2,:);
z1 = data(3,:);

x2 = data(7,:);
y2 = data(8,:);
z2 = data(9,:);

x3 = data(13,:);
y3 = data(14,:);
z3 = data(15,:);

x_mins = [min(x1),min(x2),min(x3)];
x_maxes = [max(x1),max(x2),max(x3)];
y_mins = [min(y1),min(y2),min(y3)];
y_maxes = [max(y1),max(y2),max(y3)];
z_mins = [min(z1),min(z2),min(z3)];
z_maxes = [max(z1),max(z2),max(z3)];

mins = [min(x_mins),min(y_mins),min(z_mins)];
maxes = [max(x_maxes),max(y_maxes),max(z_maxes)];

min = min(mins);
max = max(maxes);

figure;
plot3(x1,y1,z1)
hold on;
grid on;
xlabel('x')
ylabel('y')
zlabel('z')
plot3(x2,y2,z2)
plot3(x3,y3,z3)
xlim([min max])
ylim([min max])
zlim([min max])
hold off;