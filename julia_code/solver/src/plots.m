clear;
close all;

C = 1.0;		% set relativistic C
mSun = 1.0; 	% set the mass to be in solar mass units

C_CGS = 2.998e10;
G_CGS = 6.674e-8;
mSun_CGS = 1.989e33;
AU = 1.496e14; % 1 AU in cm

G = 1.0;
M = mSun_CGS;		%units of mass
L = M * (G_CGS / G) * ((C / C_CGS)^2);		%units of length
T = L * C / C_CGS;		%units of time

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

a_in = data(28,:);
e_in = data(29,:);
i_in = data(30,:);
w_in = data(31,:);
Omega_in = data(32,:);

a_in = a_in * L / AU;
a_in = a_in(find(abs(a_in) < 0.05));

a_out = data(33,:);
e_out = data(34,:);
i_out = data(35,:);
w_out = data(36,:);
Omega_out = data(37,:);

a_out = a_out * L / AU;
a_out = a_out(find(abs(a_out) < 1));

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