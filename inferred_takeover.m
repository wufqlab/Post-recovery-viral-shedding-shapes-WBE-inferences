function summary_box = black_box(pars,rel_time)

beta_E = pars(1);
beta_I = pars(2);
beta_R = pars(3);

t_E = pars(4);
t_I = pars(5);
t_R = pars(6);

%% Solve SEIR-V
% SEIR model 
% Initials 
 S0 = 1e6;
 E0 = 1;
 I0 = 1;
 R0 = 0;
 V0 = 0;
 C0 = 0;

ICs = [S0, E0, I0, R0, V0, C0, 0, 0, 0]; %the last 3 zeros correspond to the contribution from E I R

% Parameters
 lambda = 4e-7;
 k      = 1/t_E;
 delta  = 1/t_I;
 sigma  = 1/t_R;

 alpha  = 125; % 50-700
 gamma  = 0.1; % Phan et al. 2023 - range 0.03 - 0.28

 params = [lambda, k, delta, sigma,...
          beta_E, beta_I, beta_R,...
          alpha, gamma];

t_final = 360;
t_outbreak = linspace(0,t_final,t_final+1);

[T1,Y1] = ode45(@SEIRV,t_outbreak,ICs,[],params);
t = T1;

%% Examine virus contribution from each class

%%%%
%daily incidence
% tdaily  = T1(2:end);
C_cumul = Y1(:,6);
% C_daily = diff(C_cumul);

%find indexes relative to X% of final size (C_cumul(end))
final_size = C_cumul(end);

% rel_time = [0.05, 0.15, 0.25, 0.4, 0.5, 0.6, 0.75, 0.85, 0.95];
time_index = nan(size(rel_time));

for ii = 1:length(rel_time)
    
    time_index(ii) = find(C_cumul>=rel_time(ii)*final_size,1);
    
end

%daily viral shed contribution from each class and total
E_cont = diff(Y1(:,7));
I_cont = diff(Y1(:,8));
R_cont = diff(Y1(:,9));
% V_tot  = diff(Y1(:,5)); %total

%number of people in E, I, R classes.
E_count = Y1(:,2);
I_count = Y1(:,3);
R_count = Y1(:,4);

%collect
first_part_summary_box = nan(length(rel_time),7);

for jj=1:length(rel_time)

    first_part_summary_box(jj,:) = [t(time_index(jj)), E_cont(time_index(jj)), I_cont(time_index(jj)), R_cont((time_index(jj))),...
                                    E_count(time_index(jj)), I_count(time_index(jj)), R_count(time_index(jj))];

end


%% Detectability of new outbreak
measurement_uncertainty = 0; %0.25; %25% proportional error

%delay detection (I)
upper_R = R_cont*(1+measurement_uncertainty);
lower_R = R_cont*(1-measurement_uncertainty);

%%%% Estimate number of new cases needed to generate detectable increase.

VL_shed_by_I_24_hrs = (1-measurement_uncertainty)*alpha*(1-gamma)*beta_I;

required_new_cases_lb_I = max(round(measurement_uncertainty*R_cont/VL_shed_by_I_24_hrs),1); %known background
required_new_cases_ub_I = max(round(0.05*upper_R/VL_shed_by_I_24_hrs),1); %unknown background

%real-time detection (E)

VL_shed_by_E_24_hrs = (1-measurement_uncertainty)*alpha*(1-gamma)*beta_E;

required_new_cases_lb_E = max(round(measurement_uncertainty*R_cont/VL_shed_by_E_24_hrs),1); %known background
required_new_cases_ub_E = max(round(upper_R/VL_shed_by_E_24_hrs),1); %unknown background


%collect
Second_part_summary_box = nan(length(rel_time),4); %I - known/unknown background + E - known/unknown background

for jj=1:length(rel_time)

    second_part_summary_box(jj,:) = [required_new_cases_lb_I(time_index(jj)), required_new_cases_ub_I(time_index(jj)),...
                                     required_new_cases_lb_E(time_index(jj)), required_new_cases_ub_E(time_index(jj))];

end

summary_box = [first_part_summary_box, second_part_summary_box];

end






