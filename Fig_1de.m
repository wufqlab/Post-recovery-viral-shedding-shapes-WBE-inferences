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

%%
SPLIT = [78, 140];  
  
RRR = nan(5,1);
for ii=1:length(SPLIT)

%% Load data
load('data.mat')

V = cRNA2.*F2;
split = SPLIT(ii); 

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

[best_params,SSE] = run(ms,problem,500);

parameter = ["lambda";"alpha";"beta";"E(0)";"SSE"];
estimated_val = [best_params';SSE];
t = table(parameter,estimated_val);

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
RX0 = 0;

X70 = V(1);
X80 = 0;
ICs = [ICs X70 X80 RX0];

[T,Y] = ode45(@SEIRV,1:length(cRNA2(1:split)),ICs,[],best_params);


%% Plot
time = datetime(2020,10,7) + caldays(0:length(cRNA2)-1);

%% Subfigure 5
figure(ii)
    ylabel('Predicted cases','FontSize',16);
    xlabel('Reported cases','FontSize',16);
    box on; hold on;
    
    y = (diff(Y(:,9)));
    x = (newRepCases2(2:split));
    X = [ones(length(x),1) x];
    b = X\y;

    yCalc2 = X*b;%b1*x;
    scatter(x,y,20,'k','filled');
    plot(x,yCalc2,'r','LineWidth',2)

    ylim([0 inf])
    ax = gca;
    ax.YAxis.Exponent = 4;

    %calculate R2
    Rsq2 = 1 - sum((y-yCalc2).^2)/sum((y-mean(y)).^2);

    [R, P] = corrcoef(x,y); 
    
    dim = [.62 .01 .85 .3];

    correlation_R = round(R(1,2),2);
    str = strcat('R = ',{' '}, num2str(correlation_R));
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',16,'EdgeColor',[1 1 1]);
%     set(gca,'xtick',[]);
    RRR(ii) = correlation_R;

    f = gcf;

    if ii==1
        exportgraphics(f,'fig_1d.pdf','Resolution',600)
        exportgraphics(f,'fig_1d.png','Resolution',600)
    elseif ii==2
        exportgraphics(f,'fig_1e.pdf','Resolution',600)
        exportgraphics(f,'fig_1e.png','Resolution',600)
    end

end


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
    RX0 = 0;
    ICs  = [ICs X70 X80 RX0];

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
    betaR = 0.5*beta; 

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
    dy(6) = lambda*S*I;       % track cumulative cases
    dy(7) = alpha*(1-eta)*beta*I;
    dy(8) = alpha*(1-eta)*betaR*R;
    dy(9) = sigma*E;       % track recovered cases
end