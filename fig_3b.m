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
% Vals = black_box(beta_E,beta_I,beta_R,t_E,t_I,t_R);

%% set up the loop to get the value for heat map
params = [beta_E,beta_I,beta_R,t_E,t_I,t_R];
index_p1 = 6; %t_R on x-axis
index_p2 = 3; %beta_R on y-axis

N = 101; % mesh for heatmap

fold_change_p1 = [0.1,10];
fold_change_p2 = [0.1,100];

Original_p1 = params(index_p1);   %x-axis
Original_p2 = params(index_p2);   %y-axis

Range_p1    = [7, 84];
Range_p2    = log10([0.01*beta_I beta_I]);

Collection_p1 = linspace(Range_p1(1),Range_p1(2),N);
Collection_p2 = linspace(Range_p2(1),Range_p2(2),N);

% Make a cell to contain all of these
rel_time = [0.999];
num_time_interest = length(rel_time);
Collection = cell(4,num_time_interest,N,N); %4 sets of heatmaps at 9 time points each

% Now the loop
for ii=1:N

params(index_p1) = Collection_p1(ii);

    for jj=1:N

        params(index_p2) = 10^Collection_p2(jj);
        
        Vals = inferred_takeover(params,rel_time); %note that Vals contains info for 4 plots x time points
        
        % number of new cases (in I) needed with known background - slightly delay detection
        num_I_known = Vals(:,8);

        % number of new cases (in I) needed with unknown background - slightly delay detection
        num_I_unknown = Vals(:,9);

        % number of new cases (in E) needed with known background - near real time detection
        num_E_known = Vals(:,10);

        % number of new cases (in E) needed with unknown background - near real time detection
        num_E_unknown = Vals(:,11);

        % Now we collect everything in the respective position
        for kk = 1:num_time_interest

            Collection{1,kk,ii,jj} = num_I_known(kk);
            Collection{2,kk,ii,jj} = num_I_unknown(kk);
            Collection{3,kk,ii,jj} = num_E_known(kk);
            Collection{4,kk,ii,jj} = num_E_unknown(kk);

        end

    end

end

%% Plotting Heatmap
%% set 1
%rescale
min_p2 = log10(Range_p2(1)); 
max_p2 = log10(Range_p2(2));

min_p1 = Range_p1(1); 
max_p1 = Range_p1(2); 

% relative change curve
XX = Collection_p1;
Baseline_YY = Original_p1*Original_p2;
YY = Baseline_YY./XX;
YY = log10(YY);

    f = figure(1);

    Max_color = 5000;

% mm is time
for mm = 1:num_time_interest

    Matrix_holder = nan(N,N);

    for ii = 1:N
        for jj = 1:N
            Matrix_holder(ii,jj) = Collection{2,mm,ii,jj};
        end
    end

    s = pcolor(Collection_p1,Collection_p2,Matrix_holder');
    s.FaceColor = 'interp';
    set(s, 'EdgeColor', 'none')

    ylim = [Range_p2(1) Range_p2(2)];
    xlim = [min_p1 max_p1];

    yticks(log10([10^Range_p2(1), 10^Range_p2(2)/10, 10^Range_p2(2)]));
    xticks([min_p1, 28, 56, max_p1]);

    yticklabels({'1%','10%','100%'})

    if mm > 1 && mm ~= 6 
        set(gca,'ytick',[]);
    end

    if mm <= num_time_interest/2
        set(gca,'xtick',[]);
    end

    h = colorbar;
    clim([0, Max_color]); %color map min/max

    cmap = parula(25); 
    colormap(cmap)

    colorbar off;

    set(gca,'FontSize',14)
    
    hold off;

end

cbh = colorbar;
clim([0, Max_color]);
cbh.Label.String = {'Number of infectious cases needed'};
cbh.Label.FontSize = 16;
cbh.Label.Rotation = 270;
cbh.Label.Position(1) = 4.5;
cbh.Ruler.Exponent = 3;

ylabel({'Recovered shedding rate'; 'relative to infectious shedding rate'},'FontSize',16)
xlabel('Duration of shedding (days) ','FontSize',16)
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)

%% Save Figures
f2 = gcf;

exportgraphics(f2,'fig_3b.pdf','Resolution',300)
exportgraphics(f2,'fig_3b.png','Resolution',300)











