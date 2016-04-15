%% test KPP
%% 

addpath('../model_1D')

% parameters
D = .5;
N = 5000;
x = -10:.1:10;
T = 2;
dt = .01;
rhoIC = exp(-x.^2);

dx = x(2)-x(1);
2*D*dt/dx^2

% solution
rhoMacro = KPP_Macro(rhoIC,x,D,dt,T,2);
rhoMicro = KPP_Micro(rhoIC,x,D,dt,T,N);

% plot
for k = 1:length(rho(1,:))
    plot(x,rhoMicro(:,k), ...
         x,rhoMacro(:,k),'linewidth',2); 
    legend('Micro','Macro')
    title(['time t=',num2str(k*dt,'%10.2f')])
    axis([x(1) x(end) 0 1.1])
    pause(.001);
end