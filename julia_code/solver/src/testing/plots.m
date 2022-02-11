clear;
close all;

m1 = 1.4;
m2 = 1.4;

C = 1.0;		% set relativistic C
mSun = 1.0; 	% set the mass to be in solar mass units

C_CGS = 2.998e10;
G_CGS = 6.674e-8;
mSun_CGS = 1.989e33;
AU = 1.496e14; % 1 AU in cm
KM = 1e6; % 1 km in cm

G = 1.0;
M = mSun_CGS;		%units of mass
L = M * (G_CGS / G) * ((C / C_CGS)^2);		%units of length
T = L * C / C_CGS;		%units of time

tfinal = 30;        %full runtime in years

%data = importdata('PMdata.csv');
data = importdata('newtondata.csv');

t = data.data(:,1);
qx1 = data.data(:,2);
qy1 = data.data(:,3);
qz1 = data.data(:,4);
px1 = data.data(:,5);
py1 = data.data(:,6);
pz1 = data.data(:,7);

qx2 = data.data(:,8);
qy2 = data.data(:,9);
qz2 = data.data(:,10);
px2 = data.data(:,11);
py2 = data.data(:,12);
pz2 = data.data(:,13);

% Graphing results

%transform to AU units
qx1 = qx1 * L / KM;
qy1 = qy1 * L / KM;
qz1 = qz1 * L / KM;
qx2 = qx2 * L / KM;
qy2 = qy2 * L / KM;
qz2 = qz2 * L / KM;

x_mins = [min(qx1),min(qx2)];
x_maxes = [max(qx1),max(qx2)];
y_mins = [min(qy1),min(qy2)];
y_maxes = [max(qy1),max(qy2)];
z_mins = [min(qz1),min(qz2)];
z_maxes = [max(qz1),max(qz2)];

mins = [min(x_mins),min(y_mins),min(z_mins)];
maxes = [max(x_maxes),max(y_maxes),max(z_maxes)];

min_val = min(mins);
max_val = max(maxes);

figure;
plot3(qx1,qy1,qz1)
hold on;
grid on;
xlabel('x')
ylabel('y')
zlabel('z')
plot3(qx2,qy2,qz2)
xlim([min_val max_val])
ylim([min_val max_val])
zlim([min_val max_val])
hold off;

