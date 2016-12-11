% function TwoDStochasticQuick
% This function computes rhoMicro, rhoMacro for a given initial condition
% in two dimensions. It calls Micro.m and Macro.m. 

% Author: Shane Lubold (shane.lubold@asu.edu).
clc; clear all; close all;

%% Parameters
D = 1;                                % diffusion coefficients
%% Parameters micro
N = 10^2;                               % number of cells (initially)
Mass_rhoIC = 2;
seedNum = 1;
% domain
dx = .1;
dy = .1;
x = linspace(-5,5,2^5);
y = linspace(-5,5,2^5);
dt = .01;
T =   2;
t = 0:dt:T;
dr = .5;
intR = 0:dr:5;

% Initial Conditiion
mu1 = 0;
mu2 = 0;
var1 = 1;
var2 = 1;
shouldPlot = 0;
r = linspace(1.5,40,100);
VelMac = zeros(length(t),length(r));
VelMic = zeros(length(t),length(x));


[rhoMicro, Part] = Micro(x,y,D,dt,T,N,mu1,mu2,var1,var2);
[rhoDir,rho] = Macro(x,y,D,dt,T,mu1,mu2,var1,var2,true);

for l = 1:length(t)
   VelMac(l,:) = ComputeVelocity(x,y,rho(:,:,l),r); 
end
 

for q = 1:length(t)-1
% Call Micro, Macro
%------------------

% rhoMicro(:,:,q) = rhoMicro(:,:,q)/sum(sum(rhoMicro(:,:,q)));
% rho(:,:,q) = rho(:,:,q)/sum(sum(rho(:,:,q)));
Densities(q,:,:) = radial_density(x,y,rhoMicro(:,:,q), rho(:,:,q), Part,intR);  %return values of "1D" densities
WD(q) = WD_Cont(Densities(q,:,1)',Densities(q,:,2)',intR); % Compute the WD


% Compute speed of traveling waves
% ---------------------------------
%addpath('lib')
% c = .2; % threshold for computing contours.
% for j = 1:floor(T/dt + .5)
%     distanceMacro(j) = ComputeVelocity(x,y,rho(:,:,j),c);
% end
% Fit = polyfit(t(2:end),distanceMacro,1);
% figure; 
% plot(t(2:end),distanceMacro,t(2:end),polyval(Fit,t(2:end)))
% xlabel('Time'); ylabel('Velocity');
% 

   
 end
 %xval = reshape(xval,length(x)^2, length(y)^2);
 
