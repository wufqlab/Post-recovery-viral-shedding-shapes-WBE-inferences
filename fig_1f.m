close all;
clear all;

figure(1); box on; hold on;

max_range = 182;
peak = 82;

y1 = [1 1 0.99 0.98 0.95 0.94 0.94 0.95 0.96]; %To get these, run fig_1de for different cut-offs (indicated in x1)
y2 = [1 1 0.99 0.96 0.89 0.87 0.84 0.86 0.88];
x1 = [42, 60, 78, 96, 114, 120, 140, 160, 180]/max_range; %different cut-offs
x2 = x1+0.005;

xline(peak/max_range,'--r','LineWidth',2);

scatter(x1,y1,80,'LineWidth',2);
scatter(x2,y2,80,'LineWidth',2);

ylim([0.8 1])
xlim([0.2 1])

legend('','Without recovered shedding','With recovered shedding','Location','SouthWest')

xlabel('Fraction of data used')
ylabel({'Correlation coefficient'; 'viral load vs. daily incidence'})
set(gca,'FontSize',16)


f = gcf;
exportgraphics(f,'fig_1f.pdf','Resolution',600)
exportgraphics(f,'fig_1f.png','Resolution',600)