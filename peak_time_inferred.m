function summary_box = peak_time_inferred(pars)

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

%daily incidence
C_daily = diff(Y1(:,6));
[C_peak,index_C] = max(C_daily);
t_C_peak = t(index_C);

%daily viral load
V_daily = diff(Y1(:,5));
[V_peak,index_V] = max(V_daily);
t_V_peak = t(index_V);

t_dif = t_V_peak - t_C_peak;

summary_box = t_dif;

end






