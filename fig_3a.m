%% schematic
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
beta_fixed = 4.49e7;%2.4547e7;  
lb = [1E-8  100   beta_fixed  10 ];
ub = [1E-6  300   beta_fixed  5000];
p0 = [9E-8  150   beta_fixed  2000];

%% try Global Search
gs = GlobalSearch;
ms = MultiStart('Display','iter');

problem = createOptimProblem('fmincon','x0',p0,...
    'objective',@(param)obj_fun(param,tspan,V),'lb',lb,'ub',ub);

% [best_params,SSE] = run(gs,problem);

[best_params,SSE] = run(ms,problem,50);

parameter = ["lambda";"alpha";"beta";"E(0)";"SSE"];
estimated_val = [best_params';SSE];
t = table(parameter,estimated_val)

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

[T,Y] = ode45(@SEIRV,[1 250],ICs,[],best_params);


%% Plot
time = datetime(2020,10,7) + caldays(0:length(cRNA2)-1);

%% 
figure()
%     nexttile
    box on; hold on;

    plot(T,Y(:,3),'Color',[0 158 115]/255,'LineWidth',3); %I
    plot(T,Y(:,4),'Color',[204 121 167]/255,'LineWidth',3); %R

    frac_cumulative_R = Y(:,6)./Y(end,6);
    index_999 = find(frac_cumulative_R>0.995, 1);
    xline(T(index_999),'--', {'Emerging variant'}, 'LabelHorizontalAlignment', 'left', 'FontSize', 16, 'Color',[255 15 0]/255, 'LineWidth',2)

    ylabel('Count')
    xlabel('Time')
    legend('Infectious population','Recovered population','Location','southoutside', 'Orientation', 'horizontal')

    set(gca,'xtick',[]);
    set(gca,'ytick',[]);

%%
    f = gcf;
    exportgraphics(f,'fig_3a.pdf','Resolution',600)
    exportgraphics(f,'fig_3a.png','Resolution',600)

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

function dy = SEIRV(t,y,param)
    % parameters to be fit
    lambda = param(1);
    alpha = param(2);
    beta = param(3);
    betaR = 0.25*beta; 

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
    dy(6) = gamma*I; %lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
end