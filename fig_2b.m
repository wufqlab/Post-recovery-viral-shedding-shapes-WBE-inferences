%% Parameter estimation for wastewater data with temperature-adjustement
% NOTE: To remove constraint on initial data,
% replace @(param)nonlcon(param,tspan,V) with []
rng('default');

close all; 
clear all; 
clc
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex','defaultAxesFontSize',16) 
format long

bl = '#0072BD';
br = '#D95319';

%% Load data
load('data.mat')

V = cRNA2.*F2;
split = 151;

V = V(1:split);
tspan = 1:length(V);
%% Curve-fitting
options = optimoptions('fmincon','TolX',1e-12,'TolFun',1e-12,'MaxIter',50000,'MaxFunEvals',100000,'display','off');

%    lambda alpha beta  E(0)
beta_fixed = 4.49e7;
lb = [1E-8  100   beta_fixed  10 ];
ub = [1E-6  300   beta_fixed  5000];
p0 = [9E-8  150   beta_fixed  2000];

%%
figure()

for ii=[1, 4, 3] %see note in fig_2a.m
    if ii==1
        load('best_params_1_0_0.mat')
    elseif ii==2
        load('best_params_1_25_14.mat')
    elseif ii==3
        load('best_params_1_50_14.mat')
    elseif ii==4
        load('best_params_1_25_21.mat')
    elseif ii==5
        load('best_params_1_50_21.mat')
    elseif ii==6
        load('best_params_1_75_28.mat')
    end

%% Simulate with best params

alpha = best_params(2);
beta = best_params(3);

T = 15; % temp in degC; approx average of temperatures in wave 2
traveltime = 18; % hours
k = getDecay(T);

eta = 1 - exp(-k*traveltime);

% total population served by DITP
N0 = 2300000;

E0 = best_params(4);
I0 = V(1)/(alpha*beta*(1-eta));
R0 = 0;
S0 = N0 - (E0 + I0 + R0);
V0 = V(1); % use first data point 
ICs  = [S0 E0 I0 R0 V0 E0];

X70 = V(1);
X80 = 0;
ICs = [ICs X70 X80];

if ii==1
    [T,Y] = ode45(@SEIRV_1,1:length(cRNA2),ICs,[],best_params);
elseif ii==2
    [T,Y] = ode45(@SEIRV_2,1:length(cRNA2),ICs,[],best_params);
elseif ii==3
    [T,Y] = ode45(@SEIRV_3,1:length(cRNA2),ICs,[],best_params);
elseif ii==4
    [T,Y] = ode45(@SEIRV_4,1:length(cRNA2),ICs,[],best_params);
elseif ii==5
    [T,Y] = ode45(@SEIRV_5,1:length(cRNA2),ICs,[],best_params);
elseif ii==6
    [T,Y] = ode45(@SEIRV_6,1:length(cRNA2),ICs,[],best_params);
end


%% Plot
time = datetime(2020,10,7) + caldays(0:length(cRNA2)-1);

 %% Subfigure 2
%     nexttile

    if ii==1
        plot(time2(2:151),log10(cRNA2(2:151).*F2(2:151)),'.','markersize',20,'LineWidth',2,'Color',br); hold on;
        plot(time2(152:end),log10(cRNA2(152:end).*F2(152:end)),'.','markersize',20,'LineWidth',2,'Color',bl);
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0, 0.4470, 0.7410]); 
    elseif ii==2
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]); hold on;
    elseif ii==3
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0.4940  0.1840  0.5560]); hold on;
    elseif ii==4
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]); hold on;
    elseif ii==5
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0.4940, 0.1840, 0.5560]); hold on;
    elseif ii==6
        c = plot(time2(2:end),log10(diff(Y(:,5))),'LineWidth',2,'Color',[0.3000, 0.4040, 0.8060]); hold on;
    end

    ylabel('$\log_{10}$ Viral RNA Copies');

    legend('','','No recovered shedding','MRLD','HRMD','Location','SouthWest','FontSize',16)
    xlim([time2(2) time2(151)])
    xticks([time2(2) time2(94) time2(151)])
    xtickformat('MM-dd-yy')
    set(gca,'LineWidth',1.5)


end

%%
    f = gcf;
    exportgraphics(f,'fig_2b.pdf','Resolution',600)
    exportgraphics(f,'fig_2b.png','Resolution',600)

function k = getDecay(T)
    % compute temperature-adjusted decay rate of viral RNA

    tau0 = 23.76;
    Q0 = 2.5;
    T0 = 20;

    tau = tau0*Q0^(-(T - T0)/10);

    k = log(2)/tau;

end

function dy = SEIRV_1(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = 0*beta;

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);

    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/14;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_2(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .25*beta; 

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);


    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/14;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_3(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .5*beta; 

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);

    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/14;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_4(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .25*beta; 

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);

    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/21;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_5(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .5*beta;

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);

    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/21;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_6(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .5*beta; 

    dy = zeros(6,1);
    
    S = y(1);  
    E = y(2);      
    I = y(3);
    R = y(4);
%     V = y(5);

    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    sigma = 1/3;
    gamma = 1/8;
    gamma2 = 1/28;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end