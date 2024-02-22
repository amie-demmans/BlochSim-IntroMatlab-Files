%% Sample plot code to pull from
% important here is to set a large font size, apply it to labels and axis
% Make sure plot linewidth is wide enough to be visible
% Turn legend off if not needed, or make sure it is placed in a good spot

%% Part 1: line plot, label axis, change font size, add legend
fontSz = 20;

% Set up data to plot
t = [1; 2; 3];
x = linspace(0, 3*max(t), 1000);
y = exp(-x./t);

y(4,:) = y(3,:) - y(2,:);
maxContrast = max(y(4,:));
tmax = x(y(4,:) == maxContrast);

figure; hold on
for i = 1:4
    plot( x, y(i,:), 'LineWidth', 3)
end
xlabel('Time (ms)', 'FontSize',fontSz);
ylabel('M_{xy}', 'FontSize',fontSz);
ax = gca;    ax.FontSize = fontSz;
xline(tmax) % you can add reference lines with this, or yline for horizontal
legend('t1', 't2', 't3', 't3-t2','Max Contrast', 'FontSize', fontSz, 'Location','northeast');

% Set the location on the screen for the figure window to be, and how large
% it should be
set(gcf,'Position',[100 100 800 500])

%% Sample code for Scatter plots and using the subplot feature
% Add titles for each plot, and for the whole group of plots. 

% generate some data:
% Define eigenvector of diffusivity (x,y,z)
D_max = 69;
D_local = [0; D_max; 0];
D_Intermed = [2*D_max/3; D_max/3; 0];
D_spread = [D_max/3; D_max/3; D_max/3];

% The average diffusivity will be the same
DL = sum(D_local);
DI = sum(D_Intermed);
DS = sum(D_spread);

% Calculate the root mean square of this:
N = 3; % 3 points to average over
rms_L = sqrt( sum(D_local.^2));
rms_I = sqrt( sum(D_Intermed.^2));
rms_S = sqrt( sum(D_spread.^2));

% Actual plotting code
% This generates 3 separate plots in a row, in the same figure window
fontSz = 20; % font size
pointSz = 70; % point size, note the default point shape needs to be quite large

numberRows = 1;
numberColumns = 3; 
tiledlayout(numberRows,numberColumns,"TileSpacing","compact") 

% First plot
nexttile; 
x = 1;
scatter( x-0.35, DL, pointSz,'filled'); % summed diffusivity
hold on
scatter( [x-0.05,x,x+0.05], D_local, pointSz,'filled') % directionality
scatter( x+0.35, rms_L, pointSz,'filled') % RMS
title('Anisotropic', 'FontSize',fontSz)
ax = gca;    ax.FontSize = fontSz;
% since the same thing is plotted in each, only include legend once
legend( 'Summed Diffusivity', 'Each Direction', 'RMS', 'Location','east');
xlim([0.6, 1.4]); ylim([0, D_max*1.2]);
set(gca,'xticklabel',{[]})
yline(D_max);
hold off

% 2nd plot
nexttile; 
scatter( x-0.35, DI, pointSz,'filled'); % summed diffusivity
hold on
scatter( [x-0.05,x,x+0.05], D_Intermed, pointSz,'filled') % directionality
scatter( x+0.35, rms_I, pointSz,'filled') % RMS
title('Intermediate', 'FontSize',fontSz)
ax = gca;    ax.FontSize = fontSz;
xlim([0.6, 1.4]); ylim([0, D_max*1.2]);
yline(D_max);
set(gca,'xticklabel',{[]})
hold off

% 3rd plot
nexttile; 
scatter( x-0.35, DI, pointSz,'filled'); % summed diffusivity
hold on
scatter( [x-0.05,x,x+0.05], D_spread, pointSz,'filled') % directionality
scatter( x+0.35, rms_S, pointSz,'filled') % RMS
title('Isotropic', 'FontSize',fontSz)
ax = gca;    ax.FontSize = fontSz;
xlim([0.6, 1.4]); ylim([0, D_max*1.2]);
set(gca,'xticklabel',{[]})
yline(D_max);
hold off

sgtitle('RMS check', 'FontSize',fontSz)
set(gcf,'Position',[100 100 1400 500])

