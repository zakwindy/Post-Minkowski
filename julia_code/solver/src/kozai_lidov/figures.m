clear;
close all;

C = 1.0;     % set relativistic c
mSun = 1.0;  % set mass to be in solar mass units
    %values in CGS units
C_CGS = 2.998e10;
G_CGS = 6.674e-8;
mSun_CGS = 1.989e33;
AU_CGS = 1.496e13;
KM_CGS = 1e5;
    %all runs must use G = 1.0
G = 1.0;
M = mSun_CGS;        %units of mass
L = M*(G_CGS/G)*((C/C_CGS)^2);  %units of length
T = L*C/C_CGS;       %units of time

%newton files
ni1 = importdata('newtonruns/newton_i_PBB_ICL_1.csv');
ni5 = importdata('newtonruns/newton_i_PBB_ICL_5.csv');
ni9 = importdata('newtonruns/newton_i_PBB_ICL_9.csv');

nein1 = importdata('newtonruns/newton_e_in_PBB_ICL_1.csv');
nein5 = importdata('newtonruns/newton_e_in_PBB_ICL_5.csv');
nein9 = importdata('newtonruns/newton_e_in_PBB_ICL_9.csv');

%PM files
pi1 = importdata('PMruns/PM_i_PBB_ICL_1.csv');
pi5 = importdata('PMruns/PM_i_PBB_ICL_5.csv');
pi9 = importdata('PMruns/PM_i_PBB_ICL_9.csv');

pein1 = importdata('PMruns/PM_e_in_PBB_ICL_1.csv');
pein5 = importdata('PMruns/PM_e_in_PBB_ICL_5.csv');
pein9 = importdata('PMruns/PM_e_in_PBB_ICL_9.csv');

%plot the data

figure
subplot(3,2,1);
line1 = plot(ni1.data(:,1)*T/(3600),ni1.data(:,2));
hold on;
title('a_{in}= 7.48\times10^{6}km')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
line2 = plot(pi1.data(:,1)*T/(3600),pi1.data(:,2));
%legend('Newtonian Data','PM Data');
hold off;

subplot(3,2,2);
plot(nein1.data(:,1)*T/(3600),nein1.data(:,2));
hold on;
title('a_{in}= 7.48\times10^{6}km')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein1.data(:,1)*T/(3600),pein1.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;


subplot(3,2,3);
plot(ni5.data(:,1)*T/(3600),ni5.data(:,2));
hold on;
title('a_{in}= 4.675\times10^{5}km')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi5.data(:,1)*T/(3600),pi5.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(3,2,4);
plot(nein5.data(:,1)*T/(3600),nein5.data(:,2));
hold on;
title('a_{in}= 4.675\times10^{5}km')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein5.data(:,1)*T/(3600),pein5.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(3,2,5);
plot(ni9.data(:,1)*T/(3600),ni9.data(:,2));
hold on;
title('a_{in}= 2.92188\times10^{4}km')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi9.data(:,1)*T/(3600),pi9.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(3,2,6);
plot(nein9.data(:,1)*T/(3600),nein9.data(:,2));
hold on;
title('a_{in}= 2.92188\times10^{4}km')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein9.data(:,1)*T/(3600),pein9.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

hL = legend([line1,line2],{'Newtonian Data','PM Data'});

newPosition = [0.6 .95 0.2 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units',newUnits);

sgtitle('PBB ICL System')

