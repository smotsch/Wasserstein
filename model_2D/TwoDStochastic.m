% function TwoDStochastic
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

%% Parameters
D = .01;                                % diffusion coefficients
%% Parameters micro
N = 10^2;                               % number of cells (initially)
Mass_rhoIC = 2;
seedNum = 1;
% domain
dx = .2;
dy = .2;
x = -6:dx:6;
y = -5:dy:5;
% time
dt = .1;
T = 10;
t = 0:dt:T;
% Initial Conditiion
mu1 = 0;
mu2 = 0;
var1 = .1;
var2 = .02;


% Call Micro, Macro
%------------------
rhoMicro = KPP_Micro_2D(x,y,D,dt,T,N,seedNum,Mass_rhoIC,mu1,mu2,var1,var2);
rhoMacro = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,true);
    
% Compute speed of traveling waves
%---------------------------------
addpath('lib')
c = .2; % threshold for computing contours.
for j = 1:floor(T/dt + .5)
    distanceMacro(j) = ComputeVelocity(x,y,rhoMacro(:,:,j),c);
end
Fit = polyfit(t(2:end),distanceMacro,1);
figure; 
plot(t(2:end),distanceMacro,t(2:end),polyval(Fit,t(2:end)))
xlabel('Time'); ylabel('Velocity');

break
    
%Compute the WD
[WD,Diff1,Diff2] = TwoDOptimizationCode(x,y,rhoMicro(:,:,end),...
                                        rhoMacro(:,:,end))
    
  

%     subplot(2,2,1); surf(x,y,rhoMicro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, Initial');
%     subplot(2,2,2); surf(x,y,rhoMacro(:,:,1)'); colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, Initial')
%     subplot(2,2,3); surf(x,y,rhoMicro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Micro, End')
%     subplot(2,2,4); surf(x,y,rhoMacro(:,:,end)');colorbar;
%     xlabel('x'); ylabel('y'); title('Macro, End')
    
% end