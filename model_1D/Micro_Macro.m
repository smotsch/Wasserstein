% Comparison between the Micro- and Macro- dynamics of the KPP model:
%  . Micro: interacting particle systems
%         dX_t = σdB_t
%      + 'birth process' with intensity 1-ρ.
%  . Macro:
%       ∂_t ρ = DΔ_x ρ + ρ(1-ρ)
%       

% Parameters model
D = 2/3;

% Numerical parameters
dt = .01;
T  = 10;
dx = .05;
L  = 20;
method_Kpp = 2;                           % 1: explicit, 2: Cranck-Nicholson
Ninit  = 1e4;                           % number of particles initially
seedNum = 3;
% Initial condition
x = -L:dx:L;
rhoIC = normpdf(x,0,1);                               %  Initial condition 
%rhoIC = (-1<=x).*(x<=1);

% plot
shouldPlot = false;


%---------------------------------------------------------%
%---------------------------------------------------------%

%% A.1) Solve  Macro
rhoMacro = KPP_Macro(rhoIC,x,D,dt,T,method_Kpp);
%% A.2) Solve  Micro
rhoMicro = KPP_Micro(rhoIC,x,D,dt,T,Ninit,seedNum);

%% B.1) Wasserstein distance
addpath('lib')
nT = length(rhoMacro(1,:));
stock_WD = zeros(1,nT);
for k = 1:nT
    stock_WD(k) = WD_Cont(rhoMicro(:,k)',rhoMacro(:,k)',x);
end
%% B.2) Velocity traveling wave
xStar_micro = zeros(1,nT);
xStar_macro = zeros(1,nT);
nX_half = floor(length(x)/2);
for k = 1:nT
    xStar_micro(k) = interp1(rhoMicro(nX_half:end,k)',x(nX_half:end),1/2,'left');
    xStar_macro(k) = interp1(rhoMacro(nX_half:end,k)',x(nX_half:end),1/2,'left');
end

%---------------------------------------------------------%
%---------------------------------------------------------%


%-- plot
figure(1)
intT = dt*(0:(nT-1));
plot(intT,stock_WD,'r','linewidth',2)
xlabel('time')
ylabel('WD')
title('Wasserstein distance micro/macro over time')

figure(2)
if (shouldPlot)
    jumpPlot = 1;
else
    jumpPlot = nT-1;
end
for k = 1:jumpPlot:nT
    plot(x,rhoMacro(:,k),'r','linewidth',2, ...
         x,rhoMicro(:,k),'b','linewidth',2); 
    legend('Macro','Micro'); 
    title(['Time t = ',num2str(k*dt,'%10.2f')])
    grid on; 
    axis([-L L 0 1.2]);
        pause(.01);
end

figure(3)
plot(intT,xStar_macro,'r','linewidth',2, ...
     intT,xStar_micro,'b','linewidth',2)
legend('Macro','Micro','location','northwest'); 
xlabel('time t')
ylabel('position x_start')
% estimation speed
coeff = polyfit(intT(400:end),xStar_micro(400:end),1);
speedMicro = coeff(1);
coeff = polyfit(intT(400:end),xStar_macro(400:end),1);
speedMacro = coeff(1);