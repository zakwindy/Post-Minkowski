clear;
close all;

m1 = 1.4;
m2 = 1.4;
m3 = 1.4;

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

tfinal = 30;        %full runtime in years

data = importdata('PMdata.csv');
%data = importdata('newtondata.csv');

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
qx3 = data.data(:,14);
qy3 = data.data(:,15);
qz3 = data.data(:,16);
px3 = data.data(:,17);
py3 = data.data(:,18);
pz3 = data.data(:,19);

% Post process
len = length(qx1);

M_in = m1 + m2;
M_out = M_in + m3;
mu_in = m1*m2/M_in;
mu_out = M_in*m3/M_out;

x1dot = px1/m1 + px3/M_in;
y1dot = py1/m1 + py3/M_in;
z1dot = pz1/m1 + pz3/M_in;
x2dot = px2/m2 + px3/M_in;
y2dot = py2/m2 + py3/M_in;
z2dot = pz2/m2 + pz3/M_in;
x3dot = px3/m3;
y3dot = py3/m3;
z3dot = pz3/m3;

r_in_vec = zeros(3,len);
r_out_vec = zeros(3,len);

r_in_vec(1,:) = qx1-qx2;
r_in_vec(2,:) = qy1-qy2;
r_in_vec(3,:) = qz1-qz2;
COMx = (qx1*m1 + qx2*m2)/(M_in);
COMy = (qy1*m1 + qy2*m2)/(M_in);
COMz = (qz1*m1 + qz2*m2)/(M_in);
r_out_vec(1,:) = qx3 - COMx; 
r_out_vec(2,:) = qy3-COMy;
r_out_vec(3,:) = qz3-COMz;

r_in = sqrt(r_in_vec(1,:).^2 + r_in_vec(2,:).^2 + r_in_vec(3,:).^2);
r_out = sqrt(r_out_vec(1,:).^2 + r_out_vec(2,:).^2 + r_out_vec(3,:).^2);

p_in = zeros(3,len);
p_out = zeros(3,len);

p_in(1,:) = m1*x1dot;
p_in(2,:) = m1*y1dot;
p_in(3,:) = m1*z1dot;
p_out(1,:) = m3*x3dot;
p_out(2,:) = m3*y3dot;
p_out(3,:) = m3*z3dot;

rdot_in = p_in / mu_in;
rdot_out = p_out / mu_out;

v_in = sqrt(rdot_in(1,:).^2 + rdot_in(2,:).^2 + rdot_in(3,:).^2);
v_out = sqrt(rdot_out(1,:).^2 + rdot_out(2,:).^2 + rdot_out(3,:).^2);

E_in = 0.5*(v_in.^2) - G*M_in./r_in;
E_out = 0.5*(v_out.^2) - G*M_out./r_out;

a_in = -0.5*G*M_in./E_in;
a_out = -0.5*G*M_out./E_out;

rv_in_vec = zeros(3,len);
rv_out_vec = zeros(3,len);

rv_in_vec(1,:) = r_in_vec(2,:).*rdot_in(3,:)-r_in_vec(3,:).*rdot_in(2,:);
rv_in_vec(2,:) = r_in_vec(3,:).*rdot_in(1,:)-r_in_vec(1,:).*rdot_in(3,:);
rv_in_vec(3,:) = r_in_vec(1,:).*rdot_in(2,:)-r_in_vec(2,:).*rdot_in(1,:);
rv_out_vec(1,:) = r_out_vec(2,:).*rdot_out(3,:)-r_out_vec(3,:).*rdot_out(2,:);
rv_out_vec(2,:) = r_out_vec(3,:).*rdot_out(1,:)-r_out_vec(1,:).*rdot_out(3,:);
rv_out_vec(3,:) = r_out_vec(1,:).*rdot_out(2,:)-r_out_vec(2,:).*rdot_out(1,:);

rv_in = sqrt(rv_in_vec(1,:).^2 + rv_in_vec(2,:).^2 + rv_in_vec(3,:).^2);
rv_out = sqrt(rv_out_vec(1,:).^2 + rv_out_vec(2,:).^2 + rv_out_vec(3,:).^2);

i_in = acos(rv_in_vec(3,:)./rv_in);
i_out = acos(rv_out_vec(3,:)./rv_out);

e_in = sqrt(1 - rv_in.^2./(a_in*G*M_in));
e_out = sqrt(1 - rv_out.^2./(a_out*G*M_out));

nrv_in_vec = zeros(3,len);
nrv_out_vec = zeros(3,len);

nrv_in_vec(1,:) = -rv_in_vec(2,:);
nrv_in_vec(2,:) = rv_in_vec(1,:);
nrv_out_vec(1,:) = -rv_out_vec(2,:);
nrv_out_vec(2,:) = rv_out_vec(1,:);

nrv_in = sqrt(nrv_in_vec(1,:).^2 + nrv_in_vec(2,:).^2 + nrv_in_vec(3,:).^2);
nrv_out = sqrt(nrv_out_vec(1,:).^2 + nrv_out_vec(2,:).^2 + nrv_out_vec(3,:).^2);

arg1_in = nrv_in_vec(1,:)./nrv_in;
arg1_out = nrv_out_vec(1,:)./nrv_out;
Omega_in = acos(arg1_in);
Omega_out = acos(arg1_out);

arg2_in = (a_in.*(1-e_in.^2)-r_in)./(e_in.*r_in);
arg2_out = (a_out.*(1-e_out.^2)-r_out)./(e_out.*r_out);
f_in = acos(arg2_in);
f_out = acos(arg2_out);

arg3_in = (r_in_vec(1,:).*cos(Omega_in) + r_in_vec(2,:).*sin(Omega_in))./r_in;
arg3_out = (r_out_vec(1,:).*cos(Omega_out) + r_out_vec(2,:).*sin(Omega_out))./r_out;
theta_in = acos(arg3_in);
theta_out = acos(arg3_out);

w_in = theta_in - f_in;
w_out = theta_out - f_out;

% Graphing results

%transform to AU units
qx1 = qx1 * L / AU;
qy1 = qy1 * L / AU;
qz1 = qz1 * L / AU;
qx2 = qx2 * L / AU;
qy2 = qy2 * L / AU;
qz2 = qz2 * L / AU;
qx3 = qx3 * L / AU;
qy3 = qy3 * L / AU;
qz3 = qz3 * L / AU;

x_mins = [min(qx1),min(qx2),min(qx3)];
x_maxes = [max(qx1),max(qx2),max(qx3)];
y_mins = [min(qy1),min(qy2),min(qy3)];
y_maxes = [max(qy1),max(qy2),max(qy3)];
z_mins = [min(qz1),min(qz2),min(qz3)];
z_maxes = [max(qz1),max(qz2),max(qz3)];

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
plot3(qx3,qy3,qz3)
xlim([min_val max_val])
ylim([min_val max_val])
zlim([min_val max_val])
hold off;

