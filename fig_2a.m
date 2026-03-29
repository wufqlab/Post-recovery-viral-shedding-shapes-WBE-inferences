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

%% try Global Search
figure()

for ii=[1, 4, 3] %the three scenarios in main text
    %To get these, run the model fit with different assumptions
    %The first number after 1_ is the relative shedding rate for recovered, 50 means beta_R = 0.5*beta_I
    %The second number is the duration, 14 means average 14 days in R
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
    c = plot(time2(2:end),log10(diff(Y(:,6))),'LineWidth',2); 
    hc = get(c,'Color');
    hold on
    if ii==1
        p2 = plot(time2(2:end),log10(newRepCases2(2:end)),'LineWidth',2,'Color',br);
    end
    ylabel('$\log_{10}$ Daily Incidence');

    [max1,index1] = max(diff(Y(:,6))); %simulation max
    xline(time2(index1+1),'--','LineWidth',1.5,'Color',hc)
    [max2,index2] = max(newRepCases2); %case max

    if ii==1
    xline(time2(index2),'--','LineWidth',1.5,'Color',br)
    end


    ylim([2.379 4.5])
    xlim([time2(2) time2(151)])

    xticks([time2(2) time2(94) time2(151)])
    xtickformat('MM-dd-yy')
    set(gca,'LineWidth',1.5)

end

%%
    f = gcf;
    legend('','Data','Location','SouthEast','FontSize',16)

    exportgraphics(f,'fig_2a.pdf','Resolution',600)
    exportgraphics(f,'fig_2a.png','Resolution',600)


%% Functions

function err = obj_fun(param,tspan,data)
    T = 15; % temp in degC; aprrox average of temperatures in wave 2
    traveltime = 18; % hours
    k = getDecay(T);

    eta = 1 - exp(-k*traveltime);

    % total population served by DITP
    N0 = 2300000;

    E0 = param(4);
    I0 = data(1)/(param(2)*param(3)*(1-eta));
    R0 = 0;
    S0 = N0 - (E0 + I0 + R0);
    V0 = data(1);                
    ICs  = [S0 E0 I0 R0 V0 E0];

    X70 = data(1);
    X80 = 0;
    ICs  = [ICs X70 X80];

    [~,Y] = ode45(@SEIRV,tspan,ICs,[],param(1:4));

    % get daily virus
    cumVirus = Y(:,5);
    dailyVirus = diff(cumVirus);

    temp = log10(data(2:end)) - log10(abs(dailyVirus));
    adiff = rmmissing(temp);

    err = sum((adiff).^2);
end

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
    betaR = 0*beta; %10^4.67273;%10^3.966;

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
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
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
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
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
%     V = y(5);

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
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
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
%     V = y(5);

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
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
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
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end

%%
function dy = SEIRV_6(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = .75*beta;

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
    gamma2 = 1/28;
    

    dy(1) = -lambda*S*I;
    dy(2) = lambda*S*I - sigma*E;                               
    dy(3) = sigma*E - gamma*I;
    dy(4) = gamma*I - gamma2*R;
    dy(5) = alpha*(1-eta)*(beta*I+betaR*R);
    dy(6) = sigma*E; %lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end