close all;
clear all;

%% baseline parameters
omega1 = 72; %Phan et al. 2023 value: 71.97
omega2 = 4;  %Phan et al. 2023 value: 4
shedding_par = [omega1, omega2];

t_E = 3;  %average duration in E class
t_I = 8;  %average duration in I class
t_R = 30; %average duration in R class

tf = t_E + t_I + t_R;
t  = linspace(0,tf,25*tf);

f = viral_shedding_rate(t,shedding_par);

% V_E = integral(@(t) viral_shedding_rate(t,params), 0, t_E);
V_E = (omega1/2)*log((omega2^2+t_E^2)/(omega2^2));
V_I = (omega1/2)*log((omega2^2+(t_E+t_I)^2)/(omega2^2+t_E^2));
V_R = (omega1/2)*log((omega2^2+(t_E+t_I+t_R)^2)/(omega2^2+(t_E+t_I)^2));

% Find relative average viral shedding
beta_E = 10^(V_E/t_E);
beta_I = 10^(V_I/t_I);
beta_R = 10^(V_R/t_R);

%% Heatmap
% Heat map plots number of new cases required or detection vs. phi_R/phi_I
% On x axis, phi_R/phi_I increases based on t_R
% On y axis, phi_R/phi_I increases based on beta_R

%% set up the loop to get the value for heat map
params = [beta_E,beta_I,beta_R,t_E,t_I,t_R];
index_p1 = 6; %t_R on x-axis
index_p2 = 3; %beta_R on y-axis

N = 251; % mesh for heatmap (default should be 101). For test, use 5.

Original_p1 = params(index_p1);   %x-axis
Original_p2 = params(index_p2);   %y-axis

Range_p1    = [7, 84]; %Original_p1*fold_change_p1; %t_R
Range_p2    = log10([0.01*beta_I beta_I]);

Collection_p1 = linspace(Range_p1(1),Range_p1(2),N);
Collection_p2 = linspace(Range_p2(1),Range_p2(2),N);

% Make a cell to contain all of these
Collection = nan(N,N); 

% Now the loop
for ii=1:N

params(index_p1) = Collection_p1(ii);

    for jj=1:N

        params(index_p2) = 10^Collection_p2(jj);
        
        peak_time_dif = peak_time_inferred(params); %note that Vals contains info for 4 plots x time points

        Collection(ii,jj) = peak_time_dif;
       
    end

end

%% Plotting Heatmap
%% set 1 - I with known background

%rescale
min_p2 = log10(Range_p2(1));
max_p2 = log10(Range_p2(2)); 

min_p1 = Range_p1(1); 
max_p1 = Range_p1(2);

figure(1); 

    Matrix_holder = nan(N,N);

    for ii = 1:N
        for jj = 1:N
            Matrix_holder(ii,jj) = Collection(ii,jj);
        end
    end

    s = pcolor(Collection_p1,Collection_p2,Matrix_holder');
    s.FaceColor = 'interp';
    set(s, 'EdgeColor', 'none')

    ylim = [Range_p2(1) Range_p2(2)];
    xlim = [min_p1 max_p1];

    yticks(log10([10^Range_p2(1), 10^Range_p2(2)/10, 10^Range_p2(2)]));

    yticklabels({'1\%','10\%','100\%'})

    xticks([min_p1, 28, 56, max_p1]);

    h = colorbar;
    Max_color = 14;
    clim([0, Max_color]);
    cmin = 0;
    cmax = Max_color;

    colorbar off;
    set(gca,'FontSize',14)

cbh = colorbar;
clim([0, Max_color]);
cbh.Label.String = {'Peak time difference'; 'wastewater vs. daily incidence'};
cbh.Label.FontSize = 16;
cbh.Label.Rotation = 270;
cbh.Label.Position(1) = 5.2;

ylabel({'Recovered shedding rate'; 'relative to infectious shedding rate'},'FontSize',16)
xlabel('Duration of shedding (days) ','FontSize',16)
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)

%% Save Figures
f2 = gcf;
f2.Position = [100 100 700 500];
exportgraphics(f2,'fig_2c.pdf','Resolution',300)
exportgraphics(f2,'fig_2c.png','Resolution',300)












