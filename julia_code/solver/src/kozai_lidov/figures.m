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
ni1 = importdata('newtonruns/newton_i_PBSB_ICL_1.csv');
ni2 = importdata('newtonruns/newton_i_PNN_ICL_2.csv');
ni3 = importdata('newtonruns/newton_i_PNN_ICL_3.csv');
ni4 = importdata('newtonruns/newton_i_PNN_ICL_4.csv');
ni5 = importdata('newtonruns/newton_i_PBSB_ICL_4.csv');
ni6 = importdata('newtonruns/newton_i_PNN_ICL_6.csv');
ni7 = importdata('newtonruns/newton_i_PNN_ICL_7.csv');
ni8 = importdata('newtonruns/newton_i_PNN_ICL_8.csv');
ni9 = importdata('newtonruns/newton_i_PBSB_ICL_7.csv');

nein1 = importdata('newtonruns/newton_e_in_PBSB_ICL_1.csv');
nein2 = importdata('newtonruns/newton_e_in_PNN_ICL_2.csv');
nein3 = importdata('newtonruns/newton_e_in_PNN_ICL_3.csv');
nein4 = importdata('newtonruns/newton_e_in_PNN_ICL_4.csv');
nein5 = importdata('newtonruns/newton_e_in_PBSB_ICL_4.csv');
nein6 = importdata('newtonruns/newton_e_in_PNN_ICL_6.csv');
nein7 = importdata('newtonruns/newton_e_in_PNN_ICL_7.csv');
nein8 = importdata('newtonruns/newton_e_in_PNN_ICL_8.csv');
nein9 = importdata('newtonruns/newton_e_in_PBSB_ICL_7.csv');

%PM files
pi1 = importdata('PMruns/PM_i_PBSB_ICL_1.csv');
pi2 = importdata('PMruns/PM_i_PNN_ICL_2.csv');
pi3 = importdata('PMruns/PM_i_PNN_ICL_3.csv');
pi4 = importdata('PMruns/PM_i_PNN_ICL_4.csv');
pi5 = importdata('PMruns/PM_i_PBSB_ICL_4.csv');
pi6 = importdata('PMruns/PM_i_PNN_ICL_6.csv');
pi7 = importdata('PMruns/PM_i_PNN_ICL_7.csv');
pi8 = importdata('PMruns/PM_i_PNN_ICL_8.csv');
pi9 = importdata('PMruns/PM_i_PBSB_ICL_7.csv');

pein1 = importdata('PMruns/PM_e_in_PBSB_ICL_1.csv');
pein2 = importdata('PMruns/PM_e_in_PNN_ICL_2.csv');
pein3 = importdata('PMruns/PM_e_in_PNN_ICL_3.csv');
pein4 = importdata('PMruns/PM_e_in_PNN_ICL_4.csv');
pein5 = importdata('PMruns/PM_e_in_PBSB_ICL_4.csv');
pein6 = importdata('PMruns/PM_e_in_PNN_ICL_6.csv');
pein7 = importdata('PMruns/PM_e_in_PNN_ICL_7.csv');
pein8 = importdata('PMruns/PM_e_in_PNN_ICL_8.csv');
pein9 = importdata('PMruns/PM_e_in_PBSB_ICL_7.csv');

%plot the data
%{
% first create the figure
figPos = [200 200 800 500];
figure('Color', 'w', 'Position', figPos)

% next, determine how much padding you want on each side of the axes, and in
% between axes. I usually play around with these, and the figure size until
% the layout looks correct.

leftPadding = 50/figPos(3); % the space at the left of the figure
rightPadding = 25/figPos(3); % the space at the right of the figure
horizPadding = 80/figPos(3); % the space between axes (horizontally)
topPadding = 30/figPos(4); % the space at the top of the figure
bottomPadding = 50/figPos(4); % the space at the bottom of the figure
vertPadding = 120/figPos(4); % the space between axes (vertically)

% set up the grid size
nHorizAxes = 2;
nVertAxes = 9;

% figure out how big each axes should be
horizPlotSpace = 1-leftPadding-rightPadding-(nHorizAxes-1)*horizPadding;
vertPlotSpace = 1-topPadding-bottomPadding-(nVertAxes-1)*vertPadding;
width = horizPlotSpace/nHorizAxes;
height = vertPlotSpace/nVertAxes;

myAxes = zeros(nVertAxes, nHorizAxes);

% create some sample data to plot for illustrative purposes
x = linspace(0, 2*pi);
y = sin(x);

for iRow = 1:nVertAxes
    for iCol = 1:nHorizAxes
        % calculate the position
        left = leftPadding+(iCol-1)*(width+horizPadding);
        bottom = bottomPadding+(iRow-1)*(height+vertPadding);
        position = [left bottom width height];

        myAxes(iRow, iCol) = axes('Position', position);
        plot(x, y)
        xlabel('Test Label')
        ylabel('Test Label')
        title(sprintf('axes(%d, %d)', iRow, iCol))
    end
end
%}

figure
%plot the data
subplot(2,3,1);
line1 = plot(ni1.data(:,1)*T/(3600),ni1.data(:,2));
hold on;
title('a_{in}= 5\times10^{-2}AU')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
line2 = plot(pi1.data(:,1)*T/(3600),pi1.data(:,2));
%legend('Newtonian Data','PM Data');
hold off;

subplot(2,3,4);
plot(nein1.data(:,1)*T/(3600),nein1.data(:,2));
hold on;
title('a_{in}= 5\times10^{-2}AU')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein1.data(:,1)*T/(3600),pein1.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

%{
subplot(9,2,3);
plot(ni2.data(:,1)*T/(3600),ni2.data(:,2));
hold on;
title('Halved 2 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi2.data(:,1)*T/(3600),pi2.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,4);
plot(nein2.data(:,1)*T/(3600),nein2.data(:,2));
hold on;
title('Halved 2 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein2.data(:,1)*T/(3600),pein2.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,5);
plot(ni3.data(:,1)*T/(3600),ni3.data(:,2));
hold on;
title('Halved 3 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi3.data(:,1)*T/(3600),pi3.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,6);
plot(nein3.data(:,1)*T/(3600),nein3.data(:,2));
hold on;
title('Halved 3 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein3.data(:,1)*T/(3600),pein3.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,7);
plot(ni4.data(:,1)*T/(3600),ni4.data(:,2));
hold on;
title('Halved 4 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi4.data(:,1)*T/(3600),pi4.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,8);
plot(nein4.data(:,1)*T/(3600),nein4.data(:,2));
hold on;
title('Halved 4 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein4.data(:,1)*T/(3600),pein4.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;
%}
subplot(2,3,2);
plot(ni5.data(:,1)*T/(3600),ni5.data(:,2));
hold on;
title('a_{in}= 6.25\times10^{-3}AU')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi5.data(:,1)*T/(3600),pi5.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(2,3,5);
plot(nein5.data(:,1)*T/(3600),nein5.data(:,2));
hold on;
title('a_{in}= 6.25\times10^{-3}AU')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein5.data(:,1)*T/(3600),pein5.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;
%{
subplot(9,2,11);
plot(ni6.data(:,1)*T/(3600),ni6.data(:,2));
hold on;
title('Halved 6 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi6.data(:,1)*T/(3600),pi6.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,12);
plot(nein6.data(:,1)*T/(3600),nein6.data(:,2));
hold on;
title('Halved 6 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein6.data(:,1)*T/(3600),pein6.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,13);
plot(ni7.data(:,1)*T/(3600),ni7.data(:,2));
hold on;
title('Halved 7 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi7.data(:,1)*T/(3600),pi7.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,14);
plot(nein7.data(:,1)*T/(3600),nein7.data(:,2));
hold on;
title('Halved 7 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein7.data(:,1)*T/(3600),pein7.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,15);
plot(ni8.data(:,1)*T/(3600),ni8.data(:,2));
hold on;
title('Halved 8 Times')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi8.data(:,1)*T/(3600),pi8.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(9,2,16);
plot(nein8.data(:,1)*T/(3600),nein8.data(:,2));
hold on;
title('Halved 8 Times')
xlabel('Time (hours)');
ylabel('Eccentricity');
ylim([0 1])
grid on;
plot(pein8.data(:,1)*T/(3600),pein8.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;
%}
subplot(2,3,3);
plot(ni9.data(:,1)*T/(3600),ni9.data(:,2));
hold on;
title('a_{in}= 7.8125\times10^{-4}AU')
xlabel('Time (hours)');
ylabel('i (degrees)');
grid on;
plot(pi9.data(:,1)*T/(3600),pi9.data(:,2))
%legend('Newtonian Data','PM Data');
hold off;

subplot(2,3,6);
plot(nein9.data(:,1)*T/(3600),nein9.data(:,2));
hold on;
title('a_{in}= 7.8125\times10^{-4}AU')
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

sgtitle('PBSB ICL System')

